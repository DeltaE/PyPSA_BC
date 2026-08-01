from pathlib import Path

import pandas as pd

from pypsa_bc import utils

# handles the config loading centrally
from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.reporting.logger import get_logger
from pypsa_bc.studies.cascade.physics import (
    coerce_delay_hours,
    coerce_minimum_release_m3_per_hour,
)

pypsa_aparser = AttributesParser()
PLACEHOLDER_RESERVOIR_IDS = {"", "DEFAULT", "nan", "None"}


def _valid_reservoir_ids(values) -> set[str]:
    return {
        str(value).strip()
        for value in (values or [])
        if pd.notna(value) and str(value).strip() not in PLACEHOLDER_RESERVOIR_IDS
    }


def _skip_report_path() -> Path:
    return Path("logs") / "hydro_skipped_assets.csv"


def _record_skip(skipped_assets, stage, site, reason):
    if skipped_assets is None:
        return
    skipped_assets.append(
        {
            "stage": stage,
            "asset_id": site.get("asset_id", ""),
            "upper_reservoir_id": site.get("upper_reservoir_id", ""),
            "connecting_node_code": site.get("connecting_node_code", ""),
            "hydro_type": site.get("hydro_type", ""),
            "reason": str(reason),
        }
    )


def _site_var_om_cost(site):
    """Return hydro variable O&M from either legacy or current schema."""
    if "variable_om_cost_CAD_per_MWh" in site and pd.notna(site["variable_om_cost_CAD_per_MWh"]):
        return site["variable_om_cost_CAD_per_MWh"]
    if "variable_om_costs" in site and pd.notna(site["variable_om_costs"]):
        return site["variable_om_costs"]
    raise KeyError("variable_om_cost_CAD_per_MWh|variable_om_costs")


def is_terminal_stage(down_rid):
    """
    This function checks whether a cascade is a terminal cascade or not
    return (boolean): True if stage is terminal otherwise False.
    """
    terminal_list = {
        "DEFAULT",
        "Columbia (U.S.)",
        "Brilliant",
        "Lower Bonnington",
        "Upper Bonnington",
    }
    if down_rid in terminal_list:
        return True
    elif down_rid:
        return
    else:
        return False


def get_reservoir_dict(site, reservoir, inflow, res_list, bus_dict):
    """
    This function is used to return a dictionary with 7 components to be used to build
    a hydroelectric reservoirs in PyPSA.
    Unit descriptions:
    i) Flow units always in m^3 / hr for now.
    ii)

    site: hydroelectric site (contains metadata) as a pd.series.
    reservoir: reservoir for the current hydroelectric site (contains metadata) as a pd.series.
    inflow: Timeseries of the inflow for the current hydroelectric site.
    res_list: List of unique reservoirs.
    bus_dict: Dictionary for mapping node codes to their respective elc bus (min voltage assumed)
    Components:
    1) Water bus: Bus connected to discharge, spill, and connected by links to reservoir bus. (No storage on water bus)
    2) Reservoir bus: Bus directly connected to reservoir storage.
    3) Store link: Link used to move water from water bus into the reservoirs storage
    4) Release link: Link used to release water from the reservoir bus and thereby storage into the water bus for use.
    5) Reservoir store: Store used to store water for the reservoir being modelled.
    6) Inflow generator: Generator used to be source of inflow into the system.
    7) Discharge link: Link used to produces power and then release the same water to next stage.
    8) Spill link: Link use to spill water to the next stage.

    """
    # IDs
    aid = site["asset_id"]
    res_list = _valid_reservoir_ids(res_list)
    up_rid = site["upper_reservoir_id"]
    down_rid = site["lower_reservoir_id"]
    spill_down_rid = site.get("spill_lower_reservoir_id", down_rid)
    if pd.isna(spill_down_rid) or str(spill_down_rid).strip() in {"", "DEFAULT"}:
        spill_down_rid = down_rid
    travel_time_h = coerce_delay_hours(site.get("travel_time_h"))
    minimum_release = coerce_minimum_release_m3_per_hour(
        site.get("minimum_release_m3_per_s", site.get("min_water_discharge"))
    )

    # parameters
    # NOTE: Unit change should happen here if desired otherwise leave as is..
    max_storage = reservoir["max_storage"] - reservoir["min_storage"]
    # Reservoir inflows are calibrated and stored in m3/h upstream. Multiplying
    # by 3600 again would inflate only the Generator nominal capacity and obscure
    # the unit contract used by the water balance.
    max_inflow = max(inflow)

    res_dict = {}
    # 1) add water_bus
    # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
    carrier = "Water"
    class_name = "Bus"
    res_dict["water bus"] = {
        "class_name": class_name,
        "name": " ".join([aid, carrier, class_name]),
        "carrier": carrier,
    }

    # 2) add reservoir_bus
    # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
    carrier = "Water"
    class_name = "Bus"
    res_dict["reservoir bus"] = {
        "class_name": class_name,
        "name": " ".join([up_rid, carrier, class_name]),
        "carrier": carrier,
    }

    # 3) add store_link
    # name: {ASSET_ID} Store {CLASS_NAME}
    class_name = "Link"
    res_dict["store link"] = {
        "class_name": "Link",
        "name": " ".join([aid, "Store", class_name]),
        "bus0": res_dict["water bus"]["name"],
        "bus1": res_dict["reservoir bus"]["name"],
        "efficiency": 1.0,  # mass balance
        "p_nom": max_storage,  # Should be adjusted based on timestep eventually
    }

    # 4) get release_link
    # name: {ASSET_ID} Release {CLASS_NAME}
    class_name = "Link"
    res_dict["release link"] = {
        "class_name": class_name,
        "name": " ".join([aid, "Release", class_name]),
        "bus0": res_dict["reservoir bus"]["name"],
        "bus1": res_dict["water bus"]["name"],
        "efficiency": 1.0,  # mass_balance
        "p_nom": max_storage + max_inflow,
    }

    # 5) get reservoir store
    # name: {UPPER_RESERVOIR_ID} Store {CLASS_NAME}
    # UPDATE NEEDED: Need to add a starting state for the reservoirs
    class_name = "Store"
    res_dict["reservoir store"] = {
        "class_name": class_name,
        "name": " ".join([up_rid, "Reservoir", class_name]),
        "bus": res_dict["reservoir bus"]["name"],
        "e_nom": max_storage,  # units of m^3
        # "e_initial":max_storage * 0.7, # Assumed 100% filled
        "e_cyclic": True,  #
    }

    # 6) add inflow generator
    # name: {UPPER_RESERVOIR_ID} Inflow {CLASS_NAME}
    class_name = "Generator"
    res_dict["inflow generator"] = {
        "class_name": class_name,
        "name": " ".join([up_rid, "Inflow", class_name]),
        "bus": res_dict["reservoir bus"]["name"],
        "carrier": "inflow",
        "type": "Reservoir",
        "efficiency": 1.0,  # mass_balance
        "p_nom": max_inflow,  # max(inflow series)
        "p_set": inflow,
        "p_max_pu": [flow / max_inflow if flow != 0 else 0 for flow in inflow],
        "p_min_pu": [flow / max_inflow if flow != 0 else 0 for flow in inflow],
    }

    # Terminal check
    class_name = "Bus"
    carrier = "Water"
    if down_rid in res_list:  # true if not terminal
        # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
        downstream_bus = " ".join([down_rid, carrier, class_name])
    else:
        # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
        cascade_name = site["cascade_group"]
        downstream_bus = " ".join([cascade_name, carrier, class_name])
        # Add bus for the terminal reservoir
        res_dict["terminal bus"] = {
            "class_name": "Bus",
            "name": downstream_bus,
            "carrier": "Water",  # Water not water, EL
        }

        # Add spill store
        res_dict["terminal store"] = {
            "class_name": "Store",
            "name": cascade_name,
            "bus": downstream_bus,
            "e_nom": 1e15,  # The max storage needs to retain all possible water in the model horizon
            "e_cyclic": False,
        }

    if spill_down_rid == down_rid:
        spill_downstream_bus = downstream_bus
    elif spill_down_rid in res_list:
        spill_downstream_bus = " ".join([spill_down_rid, carrier, class_name])
    else:
        spill_terminal_id = str(spill_down_rid).strip()
        spill_downstream_bus = " ".join([spill_terminal_id, carrier, class_name])
        res_dict["spill terminal bus"] = {
            "class_name": "Bus",
            "name": spill_downstream_bus,
            "carrier": "Water",
        }
        res_dict["spill terminal store"] = {
            "class_name": "Store",
            "name": f"{spill_terminal_id} Terminal Store",
            "bus": spill_downstream_bus,
            "e_nom": 1e15,
            "e_cyclic": False,
        }

    # 7) get discharge link
    ## ADD code here to grab ELC bus based on the connecting_node_code
    elc_bus_name = utils.get_gen_bus(site["connecting_node_code"], bus_dict)
    ###
    q_rated = float(site["max_water_discharge"]) * 3600  # Convert from m^3/s to units of m^3 / hr
    eff_m3_to_mwhr = site["capacity"] / q_rated
    marginal_cost = _site_var_om_cost(site) * eff_m3_to_mwhr  # Needs permanent fix for costs later.
    res_dict["discharge link"] = {
        "class_name": "Link",
        "name": " ".join([aid, "Discharge Link"]),
        "bus0": res_dict["water bus"]["name"],
        "bus1": elc_bus_name,  # NEED TO CREATE connection to ELC (ELC bus name?)
        "bus2": downstream_bus,
        "marginal_cost": marginal_cost,
        "efficiency": eff_m3_to_mwhr,
        "efficiency2": 1.0,  # mass balance
        "p_nom": q_rated,  # Should be derived to ensure larger than max(inflow, spill + discharge)
        "delay2": travel_time_h,
        "cyclic_delay2": True,
        "hydro_station": aid,
        "hydro_release_role": "turbine",
        "water_output_factor": 1.0,
        "minimum_release_m3_per_hour": minimum_release,
    }
    # 8) get spill link
    if site["max_spill"] == "DEFAULT":
        spill = 9999 * 3600  # Fill a value
    else:
        spill = float(site["max_spill"]) * 3600  # Units of m^3 / hr

    if spill == 0:  # Later this should be replaced to allow 0 spill
        spill = 9999 * 3600  # Fill a value

    res_dict["spill link"] = {
        "class_name": "Link",
        "name": " ".join([aid, "Spill Link"]),
        "bus0": res_dict["water bus"]["name"],
        "bus1": spill_downstream_bus,
        "efficiency": 1.0,  # mass_balance
        "p_nom": spill,  # Should be derived to ensure larger than max(inflow, spill + discharge)
        "marginal_cost": 0.000001,  # Small number to avoid spill
        "delay": travel_time_h,
        "cyclic_delay": True,
        "hydro_station": aid,
        "hydro_release_role": "spill",
        "water_output_factor": 1.0,
        "minimum_release_m3_per_hour": minimum_release,
    }

    return res_dict


def get_ror_dict(site, ror_ts, bus_dict):
    """
    This function creates the dictionary for a single RoR facility.
    """
    # name: {ASSET_ID} RoR Generator
    # bus" {CONNECTING_NODE_CODE} ELC Bus
    name = " ".join([site["asset_id"], "RoR", "Generator"])
    elc_bus = utils.get_gen_bus(site["connecting_node_code"], bus_dict)

    # q_rated = float(site['max_water_discharge']) * 3600 # Units of m^3 / hr
    # eff_m3_to_mwhr =  site['capacity'] / q_rated
    # marginal_cost = site["variable_om_cost_USD_per_MWh"] * eff_m3_to_mwhr

    return {
        "class_name": "Generator",
        "name": name,
        "bus": elc_bus,
        "p_nom": site["capacity"],
        "type": "RoR",
        "marginal_cost": _site_var_om_cost(site),
        "p_nom_extendable": False,  # Site already built
        # "capital_cost":site[], # no applicable since built
        "p_max_pu": ror_ts.apply(lambda x: min(x / site["capacity"], 1)),
    }


def get_ror_water_dict(site, water_inflow, bus_dict, res_list=None, reservoir=None):
    """
    Creates a custom run of river asset for pypsa which is comprised of multiple assets.
    This is done to track water
    """
    ror_water_dict = {}

    aid = site["asset_id"]
    up_rid = site["upper_reservoir_id"]
    down_rid = site["lower_reservoir_id"]
    spill_down_rid = site.get("spill_lower_reservoir_id", down_rid)
    if pd.isna(spill_down_rid) or str(spill_down_rid).strip() in {"", "DEFAULT"}:
        spill_down_rid = down_rid
    travel_time_h = coerce_delay_hours(site.get("travel_time_h"))
    minimum_release = coerce_minimum_release_m3_per_hour(
        site.get("minimum_release_m3_per_s", site.get("min_water_discharge"))
    )

    # (1) Add reservoir bus where the inflow generators is attached too
    carrier = "Water"
    class_name = "Bus"
    ror_water_dict["reservoir bus"] = {
        "class_name": class_name,
        "name": " ".join([up_rid, carrier, class_name]),
        "carrier": carrier,
    }

    if reservoir is not None:
        maximum = pd.to_numeric(reservoir.get("max_storage"), errors="coerce")
        minimum = pd.to_numeric(reservoir.get("min_storage"), errors="coerce")
        usable_storage = maximum - minimum
        if pd.notna(usable_storage) and usable_storage > 0:
            ror_water_dict["reservoir store"] = {
                "class_name": "Store",
                "name": f"{up_rid} Reservoir Store",
                "bus": ror_water_dict["reservoir bus"]["name"],
                "e_nom": float(usable_storage),
                "e_cyclic": True,
            }

    # (2) Add inflow generator. The water-bus unit contract is m3/h.
    class_name = "Generator"
    max_inflow = water_inflow.max()
    ror_water_dict["inflow generator"] = {
        "class_name": class_name,
        "name": " ".join([up_rid, "Inflow", class_name]),
        "bus": ror_water_dict["reservoir bus"]["name"],
        "carrier": "inflow",
        "type": "Reservoir",
        "efficiency": 1.0,  # mass_balance
        "p_nom": max_inflow,
        "p_set": water_inflow,
        "p_max_pu": [flow / max_inflow if flow != 0 else 0 for flow in water_inflow],
        "p_min_pu": [flow / max_inflow if flow != 0 else 0 for flow in water_inflow],
    }

    # (3) add discharge link
    ## ADD code here to grab ELC bus based on the connecting_node_code
    elc_bus_name = utils.get_gen_bus(site["connecting_node_code"], bus_dict)
    class_name = "Link"

    site_capacity = site["capacity"]  # MW
    q_rated = float(site["max_water_discharge"]) * 3600  # Units of m^3 / hr

    eff_m3_to_mwhr = site_capacity / q_rated

    res_list = _valid_reservoir_ids(res_list)
    if down_rid in res_list:
        downstream_bus = " ".join([down_rid, carrier, "Bus"])
    else:
        cascade_name = site["cascade_group"]
        downstream_bus = " ".join([cascade_name, carrier, "Bus"])
        ror_water_dict["terminal bus"] = {
            "class_name": "Bus",
            "name": downstream_bus,
            "carrier": "Water",
        }
        ror_water_dict["terminal store"] = {
            "class_name": "Store",
            "name": cascade_name,
            "bus": downstream_bus,
            "e_nom": 1e15,
            "e_cyclic": False,
        }

    if spill_down_rid == down_rid:
        spill_downstream_bus = downstream_bus
    elif spill_down_rid in res_list:
        spill_downstream_bus = " ".join([spill_down_rid, carrier, "Bus"])
    else:
        spill_terminal_id = str(spill_down_rid).strip()
        spill_downstream_bus = " ".join([spill_terminal_id, carrier, "Bus"])
        ror_water_dict["spill terminal bus"] = {
            "class_name": "Bus",
            "name": spill_downstream_bus,
            "carrier": "Water",
        }
        ror_water_dict["spill terminal store"] = {
            "class_name": "Store",
            "name": f"{spill_terminal_id} Terminal Store",
            "bus": spill_downstream_bus,
            "e_nom": 1e15,
            "e_cyclic": False,
        }

    ror_water_dict["discharge link"] = {
        "class_name": class_name,
        "name": " ".join([aid, "Discharge Link"]),
        "bus0": ror_water_dict["reservoir bus"]["name"],  # res bus
        "bus1": elc_bus_name,  # elc bus
        "bus2": downstream_bus,  # downstream res
        "marginal_cost": _site_var_om_cost(site) * eff_m3_to_mwhr,
        "efficiency": eff_m3_to_mwhr,
        "efficiency2": 1.0,
        "p_nom": q_rated,
        "delay2": travel_time_h,
        "cyclic_delay2": True,
        "hydro_station": aid,
        "hydro_release_role": "turbine",
        "water_output_factor": 1.0,
        "minimum_release_m3_per_hour": minimum_release,
    }

    class_name = "Link"

    if site["max_spill"] == "DEFAULT" or pd.isna(site["max_spill"]):
        spill_m3_per_hour = 9999 * 3600
    else:
        spill_m3_per_hour = float(site["max_spill"]) * 3600
    if spill_m3_per_hour == 0:
        spill_m3_per_hour = 9999 * 3600

    ror_water_dict["spill link"] = {
        "class_name": class_name,
        "name": " ".join([aid, "Spill Link"]),
        "bus0": ror_water_dict["reservoir bus"]["name"],  # res bus
        "bus1": spill_downstream_bus,  # downstream route
        "efficiency": 1.0,
        "p_nom": spill_m3_per_hour,
        "delay": travel_time_h,
        "cyclic_delay": True,
        "hydro_station": aid,
        "hydro_release_role": "spill",
        "water_output_factor": 1.0,
        "minimum_release_m3_per_hour": minimum_release,
    }

    return ror_water_dict


def write_reservoir_dict(hydro_sites, hydro_res, res_inflows, bus_dict, cfg, skipped_assets=None):
    """
    This function creates and writes the reservoir dictionary (pickle) which contains all the necessary components
    to be read in by PyPSA to instantiate the reservoirs in PyPSA.
    """
    # Aggregation removed from here since moved to "create_hydro_assets.py"
    # subset = ["asset_id","latitude","longitude",] # Column names used to form the duplicate list
    # sum_list = ["capacity", "annual_avg_energy", "ramp_up", "ramp_down"] # List of parameters to aggregate (assume ramping applies per turbine and asset)
    temp_df = hydro_sites[
        hydro_sites["hydro_type"].str.contains("reservoir")
    ]  # Filter for reservoir hydro assets
    # agg_hydro_sites = temp_df.groupby(by="asset_id", group_keys=False).apply(lambda x:
    #                                                 hydro.merge_assets(x, subset, sum_list)) # Make sure to aggregate the assets

    # Create dictionary for json
    # Asset_ID -> components -> attributes
    res_dict = {}
    res_list = hydro_res["asset_id"].tolist()  # list of unique reservoirs modelled here
    log = get_logger("enrich_format_hydro")
    for _, site in temp_df.iterrows():
        # IDs needed to index inflow and reservoir
        aid = site["asset_id"]
        up_rid = site["upper_reservoir_id"]

        try:
            # find matching reservoirs
            reservoir = hydro_res[hydro_res["asset_id"] == up_rid].iloc[0]
            inflow = res_inflows[up_rid]
            res_dict[aid] = get_reservoir_dict(site, reservoir, inflow, res_list, bus_dict)
        except Exception as exc:
            _record_skip(skipped_assets, "reservoir", site, exc)
            log.warning(
                "Skipping hydro reservoir asset_id=%s upper_reservoir_id=%s (%s)",
                aid,
                up_rid,
                exc,
            )

    # write pickle
    out_file = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["res"]
    utils.write_pickle(res_dict, out_file)
    utils.print_update(level=3, message="Reservoir data saved to :" + out_file)


def write_ror_dict(hydro_sites, ror_series, bus_dict, cfg, skipped_assets=None):
    """
    This function writes a dictionary containing the information needed to create the components for
    existing RoR facilities in PyPSA.
    """

    ror_dict = {}
    temp_df = hydro_sites[hydro_sites["hydro_type"] == "ror"]
    log = get_logger("enrich_format_hydro")
    for _, site in temp_df.iterrows():
        aid = site["asset_id"]
        try:
            ror_ts = ror_series[aid]
            ror_dict[aid] = get_ror_dict(site, ror_ts, bus_dict)
        except Exception as exc:
            _record_skip(skipped_assets, "ror", site, exc)
            log.warning(
                "Skipping hydro RoR asset_id=%s node=%s (%s)",
                aid,
                site.get("connecting_node_code", ""),
                exc,
            )

    # write pickle
    out_file = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["ror"]
    utils.write_pickle(ror_dict, out_file)
    utils.print_update(level=3, message="RoR data saved to :" + out_file)


def write_ror_water_dict(
    hydro_sites,
    ror_series,
    reservoir_inflows,
    bus_dict,
    cfg,
    hydro_res=None,
    skipped_assets=None,
):
    """
    This function writes a dictionary containing the information needed to create the components for
    existing RoR facilities in PyPSA.
    hydro_sites:  Dataframe with all hydroelectric sites.
    ror_series: Dataframe of RoR timeseries for model horizon (hourly).
    bus_dict: Dictionary with electrical buses.
    cfg: Configuration json.
    """
    ror_dict = {}
    temp_df = hydro_sites[
        hydro_sites["hydro_type"].isin(["ror-water", "ror-water-storage"])
    ]
    res_list = _valid_reservoir_ids(hydro_sites["upper_reservoir_id"].tolist())
    log = get_logger("enrich_format_hydro")
    for _, site in temp_df.iterrows():
        aid = site["asset_id"]
        try:
            up_rid = site["upper_reservoir_id"]
            if up_rid in reservoir_inflows:
                water_inflow = reservoir_inflows[up_rid]
            else:
                # Backward-compatible fallback for stations without flow
                # observations: convert the energy profile to m3/h.
                q_rated = float(site["max_water_discharge"]) * 3600
                water_inflow = ror_series[aid] * q_rated / float(site["capacity"])
            reservoir = None
            if hydro_res is not None:
                match = hydro_res[hydro_res["asset_id"].astype(str).eq(str(up_rid))]
                if not match.empty:
                    reservoir = match.iloc[0]
            ror_dict[aid] = get_ror_water_dict(
                site,
                water_inflow,
                bus_dict,
                res_list,
                reservoir=reservoir,
            )
        except Exception as exc:
            _record_skip(skipped_assets, "ror-water", site, exc)
            log.warning(
                "Skipping hydro RoR-water asset_id=%s node=%s (%s)",
                aid,
                site.get("connecting_node_code", ""),
                exc,
            )

    # write pickle
    out_file = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["ror_water"]
    utils.write_pickle(ror_dict, out_file)
    utils.print_update(level=3, message="RoR data sevd to :" + out_file)


def main():
    """
    This script takes the hydro generation files and formats them into dataframes which will be used
    to add network components in PyPSA along with joined by their inflow or RoR power availability
    counterparts.

    This script is also responsible for aggregating the turbines.

    The inputs to this script are the following:
    1) Hydro Generation
    2) Hydro Reservoirs

    The output of this script are the following 2 dataframes:
    1) A json file of hydro assets to be modelled as RoR assets with their corresponing attributes.
    2) A json file of hydro assets to be modelled as reservoir assets with their corresponding attributes.

    Plan will be to export a json for now.. This is so there are redundant columns and avoid need of a DF per
    component.
    """
    # get the configuration file
    cfg: dict = pypsa_aparser.data_cfg
    utils.print_update(level=1, message="Formatting hydro dataset for PyPSA...")

    (start_time, end_time) = pypsa_aparser.snapshot
    utils.print_update(level=2, message=f"Snapshot extracted:{start_time},{end_time}")

    hydro_sites = pd.read_csv(cfg["output"]["create_hydro_assets"]["hydro_generation"])
    utils.print_update(level=2, message="Hydro sites loaded...")

    res_inflows = pd.read_csv(
        cfg["output"]["reservoir_inflows"]["fname"], index_col=0, parse_dates=True
    ).loc[start_time:end_time]
    utils.print_update(level=2, message="Reservoir inflows loaded...")

    ror_series = pd.read_csv(cfg["output"]["ror_ps"]["fname"], index_col=0, parse_dates=True).loc[
        start_time:end_time
    ]
    utils.print_update(level=2, message="RoR series loaded...")

    buses = pd.read_csv(cfg["output"]["base_network"] + "/buses.csv")["name"].tolist()
    utils.print_update(level=2, message="Buses loaded...")

    # (0A) Create folders if they have not been created already
    utils.ensure_path(cfg["output"]["pypsa_dict"]["folder"])

    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    utils.print_update(level=3, message="creating bus mapping dictionary...")
    bus_dict = utils.create_standard_gen_bus_map(buses)

    skipped_assets = []

    # (1) Write pickle dictionaries for the RoR facilities.
    write_ror_dict(hydro_sites, ror_series, bus_dict, cfg, skipped_assets=skipped_assets)

    # (2) Write pickle dictionaries for the reservoirs.

    hydro_res = pd.read_csv(
        cfg["output"]["create_hydro_assets"]["hydro_reservoir"]
    )  # Purely reservoir information
    utils.print_update(level=2, message="Hydro reservoirs loaded...")

    write_reservoir_dict(
        hydro_sites, hydro_res, res_inflows, bus_dict, cfg, skipped_assets=skipped_assets
    )

    # (3) Write pickle dictionaries for the RoR-Water facilities.
    write_ror_water_dict(
        hydro_sites,
        ror_series,
        res_inflows,
        bus_dict,
        cfg,
        hydro_res=hydro_res,
        skipped_assets=skipped_assets,
    )

    skip_report = _skip_report_path()
    if skipped_assets:
        skip_report.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(skipped_assets).to_csv(skip_report, index=False)
        utils.print_update(
            level=2,
            message=(
                f"Skipped {len(skipped_assets)} hydro asset(s) due to missing/invalid inputs. "
                f"Report saved to: {skip_report}"
            ),
        )
    elif skip_report.exists():
        skip_report.unlink()
        utils.print_update(level=2, message="No hydro assets skipped; removed stale skip report.")


if __name__ == "__main__":
    main()
