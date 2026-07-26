from pypsa_bc import utils
import pandas as pd
from pathlib import Path
from pypsa_bc.reporting.logger import get_logger

# handles the config loading centrally
from pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()


def _skip_report_path() -> Path:
    return Path("logs") / "tpp_skipped_assets.csv"

def get_fuel_bus(site, fuel_type, params_cfg=None):
    '''
    Returns name of the fuel bus for a specific tpp.
    '''
    if params_cfg is None:
        params_cfg = pypsa_aparser.params_cfg
    if params_cfg["workflow"]["tpp"]["gas_grid"]:
        fuel_bus = "{} {} Bus".format(site.name, fuel_type)
    else:
        fuel_bus = "Global {} Bus".format(fuel_type)

    return fuel_bus

def get_tpp_dict(site, bus_dict, tpp_gen_types, cfg):
    '''
    Creates a thermal power plant dictionary that is used to add network components of the thermal Powerplants.
    NOTE: Assume running for hourly (need to update later)
    NOTE: Fuel assumed to be in units of MMBtu
    '''
    fuel_type = tpp_gen_types[site["gen_type"]]

    fuel_bus = get_fuel_bus(site, fuel_type, pypsa_aparser.params_cfg)
    elc_bus = utils.get_gen_bus(site['connecting_node_code'], bus_dict)

    name = " ".join([site.generation_facility_code, site["gen_type"], "Link"]) # MODIFIED 2024-10-04: site.node_code -> connecting_node_code

    # Efficiency inverse of the heat rate
    eff_hr = 1. / site["heat_rate"]

    # fuel_cost + variable_cost_MWH * conversion in terms of MMBtu (var_cost * efficiency + fuel_cost)
    # NOTE: Cost for the fuel could be listed here or assigned to production method on the NG bus. 
    marginal_cost = site["average_fuel_price_CAD_per_MMBtu"] + site["variable_om_costs"] * eff_hr
    
    # Create link + store representation of generator
    # bus_name = " ".join([fuel_type, "Bus"])
    tpp_comp_dict = {"class_name":"Link",
                    "name": name,
                    "bus0": fuel_bus,
                    "bus1": elc_bus,
                    "carrier": fuel_type,
                    "efficiency":eff_hr,
                    "type" : "Thermal",
                    "ramp_limit_up":min(site["ramp_rate_percent_per_min"]*60, 1), #* site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "ramp_limit_down":min(site["ramp_rate_percent_per_min"]*60, 1), #* site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "p_nom_extendable":False,
                    "committable":pypsa_aparser.params_cfg["workflow"]["tpp"]["UC"],
                    "min_up_time":site["min_up_time_hours"],
                    "min_down_time":site["min_down_time_hours"],
                    # "ramp_limit_start_up":row["ramp_limit_start_up"], # no data atm
                    # "ramp_limit_shut_down":row["ramp_limit_shut_down"], # no data atm
                    "p_nom":site["facility_installed_capacity"] * site["heat_rate"],
                    "marginal_cost":marginal_cost,
                    # "p_min_pu":gen_params["min_plant_load"] # watch out for the forced run condition when UC is off.
                    }

    return tpp_comp_dict

def get_cogen_dict(site, gen_params, bus_dict, hist_gen, cfg):
    '''
    Creates a thermal power plant dictionary that is used to add network components of the thermal Powerplants.
    NOTE: Assume running for hourly (need to update later)
    NOTE: Fuel assumed to be in units of MMBtu
    '''
    fuel_type = "NG" # All COGEN currently NG

    # fuel_bus = get_fuel_bus(site, fuel_type, cfg)
    elc_bus = utils.get_gen_bus(site['connecting_node_code'], bus_dict)

    name = " ".join([site.name, site["gen_type"], "Generator"])

    # Efficiency inverse of the heat rate
    eff_hr = 1. / gen_params["heat_rate"]

    # fuel_cost + variable_cost_MWH * conversion in terms of MMBtu (var_cost * efficiency + fuel_cost)
    # NOTE: Cost for the fuel could be listed here or assigned to production method on the NG bus. 
    marginal_cost =  gen_params["average_fuel_price_CAD_per_MMBtu"] + gen_params["variable_om_cost_CAD_per_MWh"] * eff_hr
    

    # Find historical generation timeseries
    mask = hist_gen['Asset Short Name'] == NG_CG_CODERS_2_AESO[site.name]
    site_gen = hist_gen[mask].apply(lambda row: row['Volume'],axis=1)
    p_nom = hist_gen[mask]['Maximum Capability'].unique()[0]

    # Create link + store representation of generator
    # bus_name = " ".join([fuel_type, "Bus"])
    cogen_comp_dict = {"class_name":"Generator",
                    "name": name,
                    "bus": elc_bus,
                    # "carrier": fuel_type,
                    "efficiency":eff_hr,
                    "ramp_limit_up":min(gen_params["ramp_rate_percent_per_min"]*60, 1), #* site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "ramp_limit_down":min(gen_params["ramp_rate_percent_per_min"]*60, 1), #* site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "p_nom_extendable":False,
                    "type" : "CoGen",
                    "committable":cfg["output"]["enrich_format_tpp"]["UC"],
                    "min_up_time":gen_params["min_up_time_hours"],
                    "min_down_time":gen_params["min_down_time_hours"],
                    "p_nom":site_gen.max() + 1, # Use max value + 1
                    "marginal_cost":marginal_cost,
                    "p_set":site_gen.apply(lambda x: max(x,0.0)),
                    # "p_min_pu":site_gen.apply(lambda x: x / p_nom), 
                    # "p_max_pu":site_gen.apply(lambda x: x / p_nom)
                    }

    return cogen_comp_dict

def write_tpp_dict(tpp_assets, bus_dict, tpp_gen_types, cfg):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing vre facilities in PyPSA.
    '''
    tpp_dict = {}
    skipped_assets = []
    log = get_logger("enrich_format_tpp")

    # NEED CODE TO ADD NG BUSES.. Will need to think about this.
    # NOTE: Likley should be this eventually, since possible to have NG without tpp.
    

    # First aggregate the units
    # No reason to aggregate it seems now...
    # subset = ["connecting_node_code"]
    # sum_list = ["install_capacity_in_mw","annual_avg_energy_unit_in_gwh/y"] # Modify this later perhaps
    # tpp_agg = tpp_assets.groupby("connecting_node_code",group_keys=False).apply(lambda x:
    #                                          utils.merge_assets(x, subset, sum_list))

    # Next create the dictionaries for PyPSA instantiation
    for _,site in tpp_assets.iterrows():

        aid = site.name
        try:
            # gen_params = tpp_assets[tpp_assets["generation_type"] == site["gen_type"]].squeeze()
            tpp_dict[aid] = get_tpp_dict(site, bus_dict, tpp_gen_types, cfg)
        except KeyError as exc:
            node_code = site.get("connecting_node_code", "")
            node_suffix = "_".join(str(node_code).split("_")[1:]) if isinstance(node_code, str) else ""
            if node_suffix and node_suffix not in bus_dict:
                reason = f"missing bus mapping for node '{node_suffix}'"
            elif str(exc).strip("\"'") not in tpp_gen_types:
                reason = f"unsupported generation type '{site.get('gen_type', 'UNKNOWN')}'"
            else:
                reason = f"key error: {exc}"

            skipped_assets.append(
                {
                    "row_index": aid,
                    "asset_id": site.get("asset_id", ""),
                    "generation_facility_code": site.get("generation_facility_code", ""),
                    "connecting_node_code": node_code,
                    "gen_type": site.get("gen_type", ""),
                    "reason": reason,
                }
            )
            log.warning(
                "Skipping TPP asset row=%s asset_id=%s node=%s gen_type=%s (%s)",
                aid,
                site.get("asset_id", ""),
                node_code,
                site.get("gen_type", ""),
                reason,
            )

    # write pickle
    out_file = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["tpp"]
    utils.write_pickle(tpp_dict, out_file)
    utils.print_update(level=2,message=f"Thermal Power Plant data for PyPSA saved to: {out_file}")

    skip_report = _skip_report_path()
    if skipped_assets:
        skip_report.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(skipped_assets).to_csv(skip_report, index=False)
        utils.print_update(
            level=2,
            message=(
                f"Skipped {len(skipped_assets)} TPP asset(s) due to missing mappings. "
                f"Report saved to: {skip_report}"
            ),
        )
    elif skip_report.exists():
        skip_report.unlink()
        utils.print_update(level=2, message="No TPP assets skipped; removed stale skip report.")

def write_cogen_dict(cogen_assets, gen_generic, bus_dict, hist_gen, cfg):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing cogeneration facilities in PyPSA for Alberta
    '''
    cogen_dict = {}
    cogen_skip_list = ["AB_BKL00_GEN", "AB_KBB00_GEN"] # NOTE: Temporary as missing a historical generation in 2021 for these

    # Next create the dictionaries for PyPSA instantiation
    for _,site in cogen_assets.iterrows():

        aid = site.name
        if aid in cogen_skip_list: # skip these
            continue
        gen_params = gen_generic[gen_generic["gen_type"] == site["gen_type"]].squeeze()
        cogen_dict[aid] = get_cogen_dict(site, gen_params, bus_dict, hist_gen, cfg)

    # NOTE: Manually add SCR5 LATER
    # aid = "AB_SCR05_GEN" # Need to figure out where to locate this one 

    # write pickle
    out_file = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["cogen"]
    utils.write_pickle(cogen_dict, out_file)
    utils.print_update(level=2,message=f"Cogeneration data for PyPSA saved to: {out_file}")


def write_ff_infrastructure(tpp_gens, tpp_gen_types, cfg):
    '''
    This function is used to create a global bus and store for fossil fuels.
    Create the PyPSA components for each fuel type:
    1) Carrier
    2) Bus
    3) Store
    '''
    ff_infrastructure = {}
    ffi_list = []


    log = get_logger("enrich_format_tpp")

    for _,site in tpp_gens.iterrows():
        fuel_type = tpp_gen_types.get(site['gen_type'])
        if fuel_type is None:
            log.warning(
                "Skipping FF infrastructure for gen_type=%s (no fuel mapping)",
                site.get("gen_type", ""),
            )
            continue
        fuel_bus_name = "Global {} Bus".format(fuel_type)
        # gen_params = gen_generic[gen_generic["generation_type"] == site['gen_type']].squeeze()

        if fuel_type not in ff_infrastructure.keys():
            ff_infrastructure[fuel_type] = 1 # To denote added already
            ffi_list.append({"class_name":"Carrier",
                            "name":fuel_type,
                            "co2_emissions":site['carbon_emissions']}
                            )
                            

            ffi_list.append({"class_name":"Bus",
                            "name": fuel_bus_name,
                            "carrier":fuel_type}
                            )
        
            ffi_list.append({"class_name":'Store',
                            "name":"Global {} Store".format(fuel_type),
                            "bus":fuel_bus_name,
                            "e_nom":1e10,
                            "e_initial":1e10}
                            )
    

    out_file = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["ff_infrastructure"]
    utils.write_pickle(ffi_list, out_file)
    utils.print_update(level=2,message=f"Fossil Fuel Infrastructure for PyPSA saved to: {out_file}")

def main():
    '''
    This script creates the dictionaries needed to instantiate the thermal power plants (TPP) in PyPSA_BC.
    '''
    # Get configuration
    cfg=pypsa_aparser.data_cfg
    utils.print_update(level=1,message="Formatting thermal power dataset for PyPSA...")
    
    # gen_generic = pd.read_csv(cfg["data"]["coders"]["gen_generic"])
    # gens = pd.read_csv(cfg["data"]["coders"]["generators"])
    # hist_gen = pd.read_csv(cfg["output"]["enrich_format_tpp"]["cogen_history"],parse_dates=True, index_col=0)
    tpp_gens = pd.read_csv(cfg["output"]["create_ext_tpp_assets"]["fname"])
    utils.print_update(level=3,message=f"Thermal Power Plant Assets loaded from : {cfg['output']['create_ext_tpp_assets']['fname']}")

    buses = pd.read_csv(cfg['output']['base_network'] + "/buses.csv")['name'].tolist()
    utils.print_update(level=3,message=f"Buses loaded from : {cfg['output']['base_network'] + '/buses.csv'}")
    
    # All generation types which are thermal PP in the CODERS dataset.
    tpp_gen_types = {'NG_CT':"NG", 'NG_CC':"NG", 'gasoline_CT':"NG",
                    'oil_CT':"Oil", "coal":"Coal", 'oil_ST':"Oil",
                    'diesel_CT':"Diesel", 'coal_CCS':"Coal", 'biomass':'Biomass',
                    "biogas":"Biomass", "Coal_CCS":"Coal",
                    'NG_CG':'NG'} # Cogeneration not treated specially here


    # Determine tpp generators in BC
    # In BC all generators are "NG" however in general this is not true.
    # NOTE: For time being there will only be a single NG bus. however, in the future.
    # Buses will need to be added for each node.

    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    bus_dict = utils.create_standard_gen_bus_map(buses)
    utils.print_update(level=3,message="Bus mapping loaded...")
    
    # (1) Write pickle dictionaries for the vre assets.
    # per fuel type now to access specific costs for each
    # for region in cfg['output']['prepare_base_network']['regions']:
    #     mask = (gens["province"] == region) & (gens["gen_type"].apply(lambda x: x in tpp_gen_types ))
    #     tpp_gens = gens[mask].copy()
    #     tpp_gens.set_index('gen_node_code', inplace=True)
    write_tpp_dict(tpp_gens, bus_dict, tpp_gen_types, cfg)

    # (2)
    # Added buses, store, carrier for each unique fuel type of thermal power plants
    write_ff_infrastructure(tpp_gens, tpp_gen_types, cfg) # NOTE: Look into next 2024-09-27!!!!!!!!!!

    
if __name__ == '__main__':
    main()