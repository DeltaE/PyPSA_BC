import pypsa
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import re
from pypsa_bc import utils
from pypsa_bc.studies.cascade.constraints import (
    add_cascade_constraints,
    attach_route_release_schedules,
)
from pypsa_bc.studies.cascade.representation import aggregate_cascade_representation
from pypsa_bc.studies.cascade.uncertainty import apply_release_uncertainty_case
from pathlib import Path
import shutil
# handles the config loading centrally

import warnings
from pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()
# Suppress specific warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def _choose_solver_name() -> str:
    """Pick the first available LP/MIP solver installed in the environment."""
    solver_candidates = ["gurobi", "highs", "glpk"]
    for solver_name in solver_candidates:
        if solver_name == "gurobi":
            try:
                import gurobipy  # noqa: F401
            except Exception:
                continue
            return solver_name
        if solver_name == "highs" and shutil.which("highs"):
            return solver_name
        if solver_name == "glpk" and shutil.which("glpsol"):
            return solver_name

    raise RuntimeError(
        "No supported solver found. Install gurobi, highs, or glpk before building the network."
    )

def make_unlimited_line_capacity(network):
    utils.print_update(level=100,message="Overriding Line Capacities to Unlimited (99999)")
    for index, line in network.lines.iterrows():
        network.lines.at[index, 's_nom'] = 99999
    return network

def is_comp_in_network(comp,network):
    '''
    This function checks to see whether a component has been added already or not.
    return: boolean indicating whether a pypsa component has already been added.
    '''
    if comp['class_name'] == 'Bus':
        return comp['name'] in network.buses.index
    if comp['class_name'] == 'Link':
        return comp['name'] in network.links.index
    if comp['class_name'] == 'Store':
        return comp['name'] in network.stores.index
    if comp['class_name'] == 'Generator':
        return comp['name'] in network.generators.index
    else:
        return False
    

def add_hydro_ror_assets(network, ror_dict):
    # Add RoR assets to the mode
    for ror in ror_dict.values():
        network.add(**ror)

def add_hydro_res_assets(network, res_dict):
    # Add Reservoir assets to the model
    for res_comps in res_dict.values():
        for comp in res_comps.values():
            if not is_comp_in_network(comp,network): # Avoid duplicate error
                network.add(**comp)

def add_hydro_ror_water_assets(network, ror_water_dict):
    # Add Reservoir assets to the model
    for ror_comps in ror_water_dict.values():
        for comp in ror_comps.values():
            if not is_comp_in_network(comp,network): # Avoid duplicate error
                network.add(**comp)
                
def add_nuclear_assets(network, nuclear_dict):    
   network.add("Generator",
        name="Nuclear1",
        bus="PeaceRiver",
        p_nom=1000,  # MW
        p_nom_extendable=True,  # or True for planning models
        efficiency=0.33,  # optional, for fuel modeling (not used if marginal_cost is $/MWh)
        # fuel="uranium",  # optional, if modeling fuel flows
        marginal_cost=2.5,  # $/MWh, mostly O&M for nuclear
        capital_cost=1.3e6 * 10,  # $/MW, converted if needed
        p_min_pu=0.3,  # nuclear can't usually run below 30% output
        ramp_limit_up=0.1,  # optional ramping constraint (e.g., 10%/hr)
        ramp_limit_down=0.1,
        p_max_pu=1.0,  # or a timeseries with planned outages
    )

    

def add_wind_assets(network, wind_dict):
    # add wind farms
    for comp in wind_dict.values():
        network.add(**comp)

def add_pv_assets(network, pv_dict):
    # add pv asset
    for comp in pv_dict.values():
        network.add(**comp)

def add_tpp_assets(network, tpp_dict):
    # add tpp
    for comp in tpp_dict.values():
        network.add(**comp)

def add_cogen_assets(network, cogen_dict):
    # add cogeneration
    for comp in cogen_dict.values():
        network.add(**comp)

def add_ff_infra(network, ff_infra_list):
    # add tpp
    for comp in ff_infra_list: # actually a list
        network.add(**comp)

def get_nearest_adm_region(point,gadm_bc):
    '''
    This function determines the administrative nearest to a bus (substation).
    This is only used for substations which are not contained within a polygon.
    '''
    distances = gadm_bc.distance(point)
    min_idx = distances.idxmin()
    return gadm_bc.loc[min_idx]["NAME_2"]


def normalize_region_name(name: str) -> str:
    """Normalize administrative-region labels to the bus names used in the build."""
    return re.sub(r"[^0-9A-Za-z]+", "", name)


def repair_bus_references(network):
    """Normalize any leftover component bus references to existing network buses."""
    existing_buses = set(network.buses.index.astype(str))

    def resolve_bus_reference(bus_name):
        if not isinstance(bus_name, str):
            return bus_name
        if bus_name in existing_buses:
            return bus_name

        normalized = normalize_region_name(bus_name)
        if normalized in existing_buses:
            return normalized

        if "-" in bus_name:
            suffix = normalize_region_name(bus_name.split("-")[-1])
            if suffix in existing_buses:
                return suffix

        return bus_name

    for component_name, bus_columns in {
        "loads": ["bus"],
        "generators": ["bus"],
        "links": ["bus0", "bus1", "bus2"],
        "lines": ["bus0", "bus1"],
    }.items():
        component = getattr(network, component_name)
        for bus_column in bus_columns:
            if bus_column in component.columns:
                component[bus_column] = component[bus_column].apply(resolve_bus_reference)


def coerce_network_string_labels(network):
    """Convert Arrow-backed string labels to plain Python string/object labels."""
    component_names = [
        "buses",
        "carriers",
        "generators",
        "links",
        "loads",
        "stores",
        "lines",
        "transformers",
        "storage_units",
    ]

    for component_name in component_names:
        if not hasattr(network, component_name):
            continue

        component = getattr(network, component_name)
        if component is None or component.empty:
            continue

        component.index = component.index.astype(object)

        for column in component.columns:
            if pd.api.types.is_string_dtype(component[column]) or pd.api.types.is_object_dtype(component[column]):
                component[column] = component[column].astype(object)


def _series_numeric(df: pd.DataFrame, column_name: str) -> pd.Series:
    """Return a numeric series for an existing column or an all-NaN series when absent."""
    if column_name not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[column_name], errors="coerce")


def prepare_and_validate_line_parameters(network, min_reactance: float = 1e-4):
    """Ensure line electrical parameters are finite and physically usable before optimize()."""
    if network.lines.empty:
        raise ValueError("Network has no transmission lines after aggregation.")

    lines = network.lines

    x = _series_numeric(lines, "x")
    seg_x = _series_numeric(lines, "segment_reactance")
    reactance_ohm = _series_numeric(lines, "reactance_ohm")
    trans_seg_x = _series_numeric(lines, "Transmission_Line_Segment_Reactance")

    # Prefer existing x, then fall back to CODERS-derived reactance columns.
    candidate_x = x.copy()
    for fallback in (seg_x, reactance_ohm, trans_seg_x):
        missing_or_bad = ~np.isfinite(candidate_x) | (candidate_x <= 0)
        candidate_x.loc[missing_or_bad] = fallback.loc[missing_or_bad]

    remaining_bad = ~np.isfinite(candidate_x) | (candidate_x <= 0)
    if remaining_bad.any():
        utils.print_update(
            level=2,
            message=(
                f"Line reactance missing/invalid for {int(remaining_bad.sum())} lines; "
                f"applying floor x={min_reactance}."
            ),
        )
        candidate_x.loc[remaining_bad] = min_reactance

    lines["x"] = candidate_x.astype(float)

    # We provide explicit x/r/b/g values, so disable unresolved custom type labels.
    if "type" in lines.columns:
        lines["type"] = ""

    for col in ("r", "b", "g"):
        numeric = _series_numeric(lines, col)
        numeric = numeric.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        lines[col] = numeric.astype(float)

    if "s_nom" in lines.columns:
        s_nom = _series_numeric(lines, "s_nom")
        if (~np.isfinite(s_nom) | (s_nom <= 0)).any():
            bad_idx = lines.index[(~np.isfinite(s_nom) | (s_nom <= 0))][:10].tolist()
            raise ValueError(
                "Invalid line capacity detected (non-finite or non-positive s_nom). "
                f"Sample line IDs: {bad_idx}"
            )

    invalid_x = lines.index[(~np.isfinite(lines["x"]) | (lines["x"] <= 0))]
    if len(invalid_x) > 0:
        raise ValueError(
            "Invalid line reactance remains after preprocessing. "
            f"Sample line IDs: {invalid_x[:10].tolist()}"
        )

    utils.print_update(
        level=3,
        message=(
            "Line parameter check passed: "
            f"x in [{lines['x'].min():.6g}, {lines['x'].max():.6g}], "
            f"count={len(lines)}"
        ),
    )

def get_ab_busmap_dict(network, regional_gdf):
    # (1) Determine mapping between buses and GADM_region
    # Determine which buses line within which gadm regions
    # output: {bus:GADM_region, ...} 
    busmap_dict = {}
    # NOTE: May need to update list for ALBERTA

    for bus,row_bus in network.buses.iterrows(): # Loop over buses
        point = Point((row_bus['x'],row_bus['y']))
        row_match = regional_gdf['geometry'].apply(lambda x: x.contains(point))

        if row_match.sum() == 1:
            busmap_dict[bus] = regional_gdf[row_match]["AID"].iloc[0]
        elif row_match.sum() == 0: # Case of LB1 and JR1
            # water, discharge, release buses...
            continue
            # print("Warning: {} is not containited within any of the administrative regions!".format(bus))
        else: # Multiple matches
            print("Error: Found multiple matches for the bus named: {}".format(bus))
            exit(3)

    return busmap_dict

def get_bc_busmap_dict(network, regional_gdf):
    # (1) Determine mapping between buses and GADM_region
    # Determine which buses line within which gadm regions
    # output: {bus:GADM_region, ...} 
    busmap_dict = {}
    # NOTE: May need to update list for ALBERTA
    special_buses = ["69_LB1_DFS", "69_LB2_GSS", "69_JRI_JCT", "69_JRI_DSS",
                     '230_CMS_GSS', "138_JOR_GSS", "69_LAJ_GSS", ""] # ONYL FOR BC
    for bus,row_bus in network.buses.iterrows(): # Loop over buses
        point = Point((row_bus['x'],row_bus['y']))
        row_match = regional_gdf['geometry'].apply(lambda x: x.contains(point)).astype(bool).values

        num_matches = row_match.sum()
        if num_matches == 1:
            busmap_dict[bus] = normalize_region_name(regional_gdf[row_match]["NAME_2"].iloc[0])
        elif num_matches == 0: # Case of LB1 and JR1
            # print('ERROR no implementation for a bus which is not contained in any region!')
            # exit(3)
            if bus in special_buses:
                busmap_dict[bus] = normalize_region_name(get_nearest_adm_region(point, regional_gdf))
            else:
                # water, discharge, release buses...
                continue
            # print("Warning: {} is not containited within any of the administrative regions!".format(bus))
        else: # Multiple matches
            print("Error: Found multiple matches for the bus named: {}".format(bus))
            exit(3)

    return busmap_dict


def get_single_region_busmap_dict(network, gadm_canada):
    '''
    (1) Determine mapping between buses and GADM_region
    Determine which buses line within which gadm regions
    output: {bus:GADM_region, ...} 
    '''
    
    busmap_dict = {}
    special_buses = ["69_LB1_DFS", "69_LB2_GSS", "69_JRI_JCT", "69_JRI_DSS"]
    for bus,row_bus in network.buses.iterrows(): # Loop over buses
        point = Point((row_bus['x'], row_bus['y']))
        row_match = gadm_canada['geometry'].apply(lambda x: x.contains(point))

        if row_match.sum() == 1:
            busmap_dict[bus] = "BC"
        elif row_match.sum() == 0: # Case of LB1 and JR1
            if bus in special_buses:
                busmap_dict[bus] = "BC"
            else:
                # water, discharge, release buses...
                continue
            # print("Warning: {} is not containited within any of the administrative regions!".format(bus))
        else: # Multiple matches
            print("Error: Found multiple matches for the bus named: {}".format(bus))
            exit(3)



        # for idx, row_gadm in gadm_bc.iterrows():
        #     if row_gadm.geometry.contains(point): # Check to see which GADM region its contained within
        #         # NOTE: Need to update this to make sure LB1 and JR1 are both assigned to regions
        #         busmap_dict[bus] = row_gadm["NAME_2"]
        #         break
    
    # for bus,row in network.buses.iterrows(): # Loop over buses
    #     if bus not in busmap_dict.keys():
    #         print("Warning: {} is not containited within any of the administrative regions!".format(bus))

    return busmap_dict

def create_adm_buses(network, gadm_bc, busmap_dict):
    '''
    Created buses for the administrative regions.
    NOTE: This will need to be updated to make better use of the voltage aggregation.
    '''
    # (2) Create a new bus for GADM regions which have been mapped
    for gadm_region in busmap_dict.values(): 
        if gadm_region not in network.buses.index.to_list():
            # Could also use a dictionary to provide these values
            region_mask = gadm_bc["NAME_2"].apply(normalize_region_name) == gadm_region
            network.add(class_name="Bus",
                        name=gadm_region,
                        x=gadm_bc[region_mask].geometry.centroid.x.iloc[0],
                        y=gadm_bc[region_mask].geometry.centroid.y.iloc[0],
                        v_nom = 300 # voltage assumed
                        )
        else:
            pass

def create_ab_pr_buses(network, regional_gdf, busmap_dict):
    '''
    Created buses for the AB planning regions. (42 of them total)
    NOTE: This will need to be updated to make better use of the voltage aggregation.
    '''
    # (2) Create a new bus for GADM regions which have been mapped
    for planning_region in busmap_dict.values():
        if planning_region not in network.buses.index.to_list():
            # Could also use a dictionary to provide these values
            network.add(class_name="Bus",
                        name=planning_region,
                        x=regional_gdf[regional_gdf["AID"] == planning_region].geometry.centroid.x.iloc[0],
                        y=regional_gdf[regional_gdf["AID"] == planning_region].geometry.centroid.y.iloc[0],
                        v_nom = 300 # voltage assumed NOTE: Needs updating... later
                        )
        else:
            pass

def create_trade_buses(network):
    '''
    Add buses to represent the US and AB.
    '''
    # busmap_dict: old_bus -> trade_bus
    # 
    network.add(class_name="Bus",
                name="US Trade",
                x=-122.873948,	
                y=48.970000,
                v_nom = 300 # voltage assumed
                )
    
    network.add(class_name="Bus",
                name="AB Trade",
                x=-114.08,
                y=49.500543,
                v_nom = 300 # voltage assumed
                )


def replace_bus_refs(component, col, busmap_dict):
    '''
    Replaces bus reference in PyPSA components.
    Any component with a matching component is mapped to new component name.
    Example: bus = 69_JRI_DSS -> ADMINISTRAIVE_REGION
    '''
    mask = component[col].isin(busmap_dict.keys())
    component.loc[mask,col] = component.loc[mask,col].map(busmap_dict)

def remove_old_components(network, busmap_dict):
    '''
    Removes the old components from the detailed network
    before it was clustered.
    '''
    # (4) remove old buses,  remove lines inside a new node

    # (i) Remove interior lines to new nodes/buses
    lines_to_remove = []
    for line_name,row in network.lines.iterrows():
        if row['bus0'] == row['bus1']:
            lines_to_remove.append(line_name)

    for line_name in lines_to_remove:
        network.remove(class_name='Line',name=line_name)

    # (ii) Remove old buses
    for bus_name in busmap_dict.keys():
        if bus_name in network.buses.index:
            network.remove(class_name='Bus',name=bus_name)

    # (iii) remove transformers
    trans_to_remove = []
    for trans_name,row in network.transformers.iterrows():
        trans_to_remove.append(trans_name)

    for trans_name in trans_to_remove:
        network.remove(class_name='Transformer',name=trans_name)

    # NOTE: Custom revmoval of Site C
    # Will need to be updated
    # stc_link = 'BC_STC_GSS Discharge Link'
    # if stc_link in network.links.index:
    #     network.remove(class_name='Link', name=stc_link)

def add_trade(network, cfg):
    '''
    This function adds load to the network.
    '''
    year = cfg['trade']['year']
    start_time = cfg['scope']['temporal']['start'] 
    end_time = cfg['scope']['temporal']['end']

    temp = pd.read_excel(cfg['trade']['path'],skiprows=1)
    temp = temp[["US Timelines", "AB Timelines"]]
    temp['TIME'] = pd.date_range(start=str(year)+'-01-01 00:00:00', end=str(year)+'-12-31 23:00:00', freq='h')
    trade = temp.copy()
    trade = trade.set_index('TIME')
    trade = trade.rename(columns={'US Tielines': 'US Trade', 'AB Tielines': 'AB Trade'})
    trade = trade.loc[start_time:end_time]
    for col in trade.columns: 
        load_ts= trade[col]
        network.add("Load", "{} ELC Load".format(col), bus=col, p_set=load_ts)


# def add_pypsa_dict(network, comp_list):
#     for comp_dict in comp_list:
#         if comp_dict['bus'] == 'CentralCoast':
#             continue
#         if comp_dict['bus'] == 'Stikine':
#             continue
#         if comp_dict['bus'] == "NorthernRockies":
#             continue
#         network.add(**comp_dict)

def add_vre_expansion_sites(network, 
                            sites, 
                            ts, 
                            vre_type,
                            capacity_choice):
    """
    Temporary function to add pypsa components of VRE expansion sites.
    
    Args:
        network: pypsa network
        sites(pd.Dataframe) ; Sites data
        ts (pd.Dataframe) : Resource Timeseries (hourly)
        capacity_choice (str): 'investment' or 'full_potential'
        vre_type (str) : 'Wind' or 'Solar'
    """

    CAD_2_USD = 1.3 # Same ratio used by the CODERS data NOTE: To be updated at a later date!
    for idx,row in sites.iterrows():
        name_id = idx # ID of the new sites
        
        # # Add to network.
        # if vre_type == "Wind":
        #     marginal_cost = 0.001
        # elif vre_type == "PV":
        #     marginal_cost = 0.001
        capacity_choice_mapping = {'investment': 'invested_capacity_MW',
                                    'full_potential':'potential_capacity'}
        network.add(
            class_name = "Generator",
            name = "New {} {}".format(vre_type, name_id),
            bus = ''.join(row['Region'].split(' ')), # Currently removes space before creating uniqut bus name
            p_max_pu = ts[name_id],
            p_nom = row[capacity_choice_mapping[capacity_choice]],
            # p_nom = row['clews_investment_MW'],
            marginal_cost = row['vom'], # NOTE: Needs to synchronized
            capital_cost = row['capex'] * CAD_2_USD * 1e6, # NOTE: Currently converting M$ USD to CAD $
            p_nom_extendable = False,
            p_nom_max = row["potential_capacity"]
            )

def add_nuc_expansion_sites(network, 
                            sites,
                            capacity_choice,
                            bus_fix='PeaceRiver'):
    """
    Temporary function to add pypsa components of VRE expansion sites.
    
    Args:
        network: pypsa network
        sites(pd.Dataframe) ; Sites data
        bus_fix (str): Bus name to be used for the new sites. Can be changed to any supported bus name in the network.
        capacity_choice (str): 'investment' or 'full_potential'
    """

    CAD_2_USD = 1.3 # Same ratio used by the CODERS data NOTE: To be updated at a later date!
    for idx,row in sites.iterrows():
        name_id = idx # ID of the new sites
        
        capacity_choice_mapping = {'investment': 'invested_capacity_MW',
                                    'full_potential':'potential_capacity'}
        network.add(
            class_name = "Generator",
            name = "New {} {}".format('Nuclear', name_id),
            type="Nuclear Power Plant",
            # bus = ''.join(row['Region'].split(' ')), # Currently removes space before creating uniqut bus name
            bus = bus_fix, # Testing for Different Regions for now
            p_max_pu = 0.9,
            p_min_pu=0.3,  # nuclear can't usually run below 30% output
            p_nom = row[capacity_choice_mapping[capacity_choice]],
            ramp_limit_up=0.1,  # optional ramping constraint (e.g., 10%/hr)
            ramp_limit_down=0.1,
            marginal_cost = row['vom'], # NOTE: Needs to synchronized
            capital_cost = row['capex'] * CAD_2_USD * 1e6, # NOTE: Currently converting M$ USD to CAD $
            p_nom_extendable = False,
            # p_nom_max = row["potential_capacity"]
            )

def add_vre_committed_sites(network, 
                            sites, 
                            ts, 
                            vre_type):
    """
    Temporary function to add pypsa components of VRE expansion sites.
    
    Args:
        network: pypsa network
        sites(pd.Dataframe) ; Sites data
        ts (pd.Dataframe) : Resource Timeseries (hourly)
        capacity_choice (str): 'investment' or 'full_potential'
        vre_type (str) : 'Wind' or 'Solar'
    """
    CAD_2_USD = 1.3 # Same ratio used by the CODERS data NOTE: To be updated at a later date!
    
    for idx,row in sites.iterrows():
        name_id = idx # ID of the new sites
        
        # # Add to network.
        # if vre_type == "Wind":
        #     marginal_cost = 0.001
        # elif vre_type == "PV":
        #     marginal_cost = 0.001
 
        network.add(
            class_name = "Generator",
            name = "CFP24 {} {}".format(vre_type, name_id),
            bus = ''.join(row['Region'].split(' ')), # Currently removes space before creating uniqut bus name
            p_max_pu = ts[name_id],
            p_nom = row['potential_capacity'],
            # p_nom = row['clews_investment_MW'],
            marginal_cost = row['vom'], # NOTE: Needs to synchronized
            capital_cost = row['capex'] * CAD_2_USD * 1e6, # NOTE: Currently converting M$ USD to CAD $
            p_nom_extendable = False,
            p_nom_max = row["potential_capacity"]
            ) 
def aggregate_lines(n):
    '''
    This function aggregates the transmission lines between administrative zones of BC into a single lines.
    For now the capacity of the lines is simply aggregated.
    
    NOTE: Aggregation needs to be modified eventually to adjust other parameters affecting the admittance.
    It also removes the old lines which were aggregated.
    
    n: PyPSA Network object
    
    '''
    temp_dict = {}

    # 1) Store aggregate results
    for _,row in n.lines.iterrows():
        name = "-".join(sorted([row['bus0'], row['bus1']]))

        if name not in temp_dict.keys():
            row_dict = row.to_dict().copy()
            row_dict["class_name"] = "Line"
            row_dict["name"] = name
            temp_dict[name] = row_dict

        else:
            row_dict = row.to_dict().copy()
            temp_dict[name]["s_nom"] +=  row_dict['s_nom']

    # 2) Remove all old lines
    for line_name in n.lines.index.to_list():
        n.remove(class_name='Line',name=line_name)
        utils.print_update(level=4,message=f"Removed line: {line_name}")

    # 3) Add new lines
    for line in temp_dict.values():
        n.add(**line)
        # utils.print_update(level=4,message=f"Added line: {line}")

def fix_vre(n, hist_gen, ab_gen_map, vre_sel):
    '''
    Function to fix the wind generation to historical levels in 2021 AB.
    vre_sel: Is either 'wind' or 'solar' 
    # Strange AB_WHL01_GEN is missing...
    '''
    for gen_code, hist_gen_code in ab_gen_map[vre_sel].items():


        gen_mask = n.generators.index.str.contains(gen_code)
        generator = n.generators.loc[gen_mask,:]
    
        if generator.shape[0] == 0:
            print('{} missing!'.format(gen_code))
            continue
        
        name = generator.index[0]

        if gen_code == "": # Case where an asset needs to be added

            mask = hist_gen['Asset Short Name'] == hist_gen_code
            p_nom_max = hist_gen.loc[mask,:]['Volume'].max() + 0.05
            p_max_pu = pd.Series([1.0]*n.snapshots.shape[0], index=n.snapshots) # Give it ability to always

            name_id = hist_gen.loc[mask,:]['Planning Asset Short Name'][0]
            bus_id = "AREA" + str(hist_gen.loc[mask,:]['Planning Area'][0])

            n.add(
            class_name = "Generator",
            name = "{} {} Generator".format(name_id, vre_sel.capitalize()),
            bus = bus_id,
            p_max_pu = p_max_pu,
            p_set = hist_gen.loc[mask,:]['Volume'],
            p_nom = p_nom_max,
            marginal_cost = 0.0001, # NOTE: Needs to synchronized
            )

            continue

        if hist_gen_code == "":
            # Remove from network
            # (ii) Remove old buses
            n.remove(class_name='Generator',name=name)

        else:
            # Initialize
            p_max_pu = pd.Series([1.0]*n.snapshots.shape[0], index=n.snapshots) # Give it ability to always produce for time being
            p_set = pd.Series([0.0]*n.snapshots.shape[0], index=n.snapshots) # Modify to len(n.snapshots)

            if type(hist_gen_code) == list: # multiple historical generations for a single assset  # noqa: E721
                print('Doing the X-to-1 asset {}'.format(hist_gen_code))
                p_nom_max = 0.05
                
                for sub_gen_code in hist_gen_code:
                    mask = hist_gen['Asset Short Name'] == sub_gen_code
                    p_nom_max += hist_gen.loc[mask,:]['Volume'].max()
                    p_set += hist_gen.loc[mask,:]['Volume']
            else:
                mask = hist_gen['Asset Short Name'] == hist_gen_code
                p_nom_max = hist_gen.loc[mask,:]['Volume'].max() + 0.05
                # p_max_pu = pd.Series([1.0]*n.snapshots.shape[0],index=n.snapshots) # Give it ability to always produce for time being
                p_set += hist_gen.loc[mask,:]['Volume']
                # temp_series = pd.Series([0.0]*n.snapshots.shape[0],index=n.snapshots)
                # p_set = temp_series.add(p_set_temp, fill_value=0.0)

            
            # Reset the following
            # 1) p_max_pu
            # 2) p_set
            # 3) p_min_pu
            
            n.generators.loc[gen_mask,'p_nom'] = p_nom_max
            n.generators_t.p_set[name] = p_set
            n.generators_t.p_max_pu[name] = p_max_pu



def fix_ng_gen(n, ab_gen_map, ng_prices, ng_dict, ng_sel):
    '''
    Function to fix the wind generation to historical levels in 2021 AB.
    ng_sel: Is either 'ngcc' or 'ngct' 
    # Strange AB_WHL01_GEN is missing...
    '''
    for gen_code,hist_gen_code in ab_gen_map[ng_sel].items():


        gen_mask = n.links.index.str.contains(gen_code)
        generator = n.links.loc[gen_mask,:]
    
        if generator.shape[0] == 0:
            print('{} missing!'.format(gen_code))
            continue
        
        name = generator.index[0]


        if hist_gen_code == '':
            # Remove from network
            # (ii) Remove old buses
            n.remove(class_name='Link',name=name)
        else:
            # 
            
            MMbtu_per_GJ = 0.947817
            temp_series = pd.Series([0.0]*8760,index=n.snapshots)
            temp_series['2022-01-01':'2022-01-31'] = ng_prices.loc['2022-01-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-02-01':'2022-02-28'] = ng_prices.loc['2022-02-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-03-01':'2022-03-31'] = ng_prices.loc['2022-03-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-04-01':'2022-04-30'] = ng_prices.loc['2022-04-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-05-01':'2022-05-31'] = ng_prices.loc['2022-05-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-06-01':'2022-06-30'] = ng_prices.loc['2022-06-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-07-01':'2022-07-31'] = ng_prices.loc['2022-07-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-08-01':'2022-08-31'] = ng_prices.loc['2022-08-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-09-01':'2022-09-30'] = ng_prices.loc['2022-09-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-10-01':'2022-10-31'] = ng_prices.loc['2022-10-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-11-01':'2022-11-30'] = ng_prices.loc['2022-11-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-12-01':'2022-12-31'] = ng_prices.loc['2022-12-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series += ng_dict[ng_sel]["variable_om_cost_CAD_per_MWh"] * (1 / ng_dict[ng_sel]["heat_rate_MMBtu_per_MWh"])
            
            n.links_t.marginal_cost[name] = temp_series

def fix_cogen(n, ab_gen_map, ng_prices, ng_dict):
    '''
    Function to fix the wind generation to historical levels in 2021 AB.
    ng_sel: Is either 'ngcc' or 'ngct' 
    # Strange AB_WHL01_GEN is missing...
    '''
    for gen_code,hist_gen_code in ab_gen_map['cogen'].items():


        gen_mask = n.generators.index.str.contains(gen_code)
        generator = n.generators.loc[gen_mask,:]
    
        if generator.shape[0] == 0:
            print('{} missing!'.format(gen_code))
            continue
        
        name = generator.index[0]


        if hist_gen_code == '':
            # Remove from network
            # (ii) Remove old buses
            n.remove(class_name='Generator',name=name)
        else:
            # 
            # print(name)
            MMbtu_per_GJ = 0.947817
            temp_series = pd.Series([0.0]*8760,index=n.snapshots)
            temp_series['2022-01-01':'2022-01-31'] = ng_prices.loc['2022-01-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-02-01':'2022-02-28'] = ng_prices.loc['2022-02-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-03-01':'2022-03-31'] = ng_prices.loc['2022-03-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-04-01':'2022-04-30'] = ng_prices.loc['2022-04-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-05-01':'2022-05-31'] = ng_prices.loc['2022-05-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-06-01':'2022-06-30'] = ng_prices.loc['2022-06-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-07-01':'2022-07-31'] = ng_prices.loc['2022-07-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-08-01':'2022-08-31'] = ng_prices.loc['2022-08-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-09-01':'2022-09-30'] = ng_prices.loc['2022-09-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-10-01':'2022-10-31'] = ng_prices.loc['2022-10-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-11-01':'2022-11-30'] = ng_prices.loc['2022-11-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series['2022-12-01':'2022-12-31'] = ng_prices.loc['2022-12-01',"NG Price (CAD/GJ)"] / MMbtu_per_GJ
            temp_series += ng_dict['cogen']["variable_om_cost_CAD_per_MWh"] * (1 / ng_dict['cogen']["heat_rate_MMBtu_per_MWh"])
            
            n.generators_t.marginal_cost[name] = temp_series


def remove_coal_gen(n):
    '''
    Custom function to remove coal generators for 2022 model run.
    '''
    coal_keepers = ['AB_GNC01_GEN Coal Link', 'AB_GNC02_GEN Coal Link', 'AB_GNC03_GEN Coal Link']
    coal_mask = n.links.index.str.contains('Coal')
    coal_names = n.links[coal_mask].index.to_list()
    for cn in coal_names:
        if cn not in coal_keepers:
            n.remove(class_name='Link', name=cn)


def add_reserves(n):
    '''
    Adds reserve formulation as seen in Bethany Frew's work.
    '''
    # Adding links
    cc_mask = n.links.index.str.contains('NG_CC') 
    ct_mask = n.links.index.str.contains('NG_CT')

    # accessing
    # var1 = n.model['Link-p']['2022-01-01 00:00:00', 'AB_AFG01_GEN Biomass Link']

    # Adding spinning reserve variable
    spin_res_var = n.model.add_variables(0, np.inf, coords=[n.snapshots, n.links.index[cc_mask]], name="Link-spin-r")
    nonspin_res_var = n.model.add_variables(0, np.inf, coords=[n.snapshots, n.links.index[ct_mask]], name="Link-nonspin-r")


    # params
    cc_eff = n.links[cc_mask]['efficiency'].unique()[0]
    ct_eff = n.links[ct_mask]['efficiency'].unique()[0]

    reg_reserve = n.loads_t.p_set.sum(axis=1).apply(lambda x: 0.015*x)
    contingency = n.loads_t.p_set.sum(axis=1).apply(lambda x: max(0.015*x,466)) # Most severre single contingency (AESO doc)
    for_reserve = 50

    # Spinning reserve constraint
    lhs = spin_res_var.loc[:,cc_mask].sum('Link') * cc_eff
    rhs = reg_reserve + 0.5 * contingency  + (1/3) * for_reserve
    spin_con = n.model.add_constraints(lhs >= rhs, name='Spinning-r')


    # Non-spinning reserve
    lhs = nonspin_res_var.loc[:,ct_mask].sum('Link') * ct_eff
    rhs = 0.5 * contingency  + (2/3) * for_reserve
    nonspin_con = n.model.add_constraints(lhs >= rhs, name='Nonspinning-r')

    # Capacity limitation (spinning)
    indices = n.links[cc_mask | ct_mask].index
    values = [n.links.loc[ind,'p_nom'] for ind in indices]
    p_nom_series = pd.Series(values,index=indices)


    lhs = spin_res_var.loc[:,cc_mask] + nonspin_res_var.loc[:,ct_mask] + n.model.variables['Link-p'].loc[:,cc_mask | ct_mask]
    rhs = p_nom_series
    capacity_lim = n.model.add_constraints(lhs <= rhs, name='Capacity-r')

def main(copperplate:bool=False,
         capacity_choice:str='investment',
         tx_line_infinity:bool=False,
         year:int=2021,
         solved_network_save_to:Path | None=None,
         include_vre_investments:bool=True,
         reservoir_representation:str="A_full_cascade",
         water_policy:str="evidence_based",
         release_multiplier:float=1.0,
         release_schedule_path:Path | None=None,
         release_uncertainty_protocol_path:Path | None=None,
         cascade_topology_path:Path | None=None,
         build_report_save_to:Path | None=None,
         solve_network:bool=True):
    '''
    This script is used to build the model.(Currently, designed to build the existing electricity system in BC. (with site-c))
    The scripts takes in the following data:
    Args:
        copperplate (bool) : 'True ' or 'False
        capacity_choice (str): 'investment' or 'full_potential'
        tx_line_infinity (Bool) : 'True' or 'False'
        year (int): 2021 to 2050
    
    1) Network structure
    2) Hydro assets
    3) Hydro timeseries data
    4) Wind assets
    5) Wind timeseries data
    6) Solar assets
    7) Solar timeseries data

    The script will save a NetCDF file of the instantiated network. 

    '''
    utils.print_update(level=1,message=" PyPSA model building initiated...")
    utils.print_update(level=100,message="Disclaimer: This model supports upto 28 Regional Districts (administrative regions) as nodes. The detailed loads (if provided) will be aggregated to these regional nodes.")
    # (0) Load config file
    
    # Print initialization summary
    utils.print_update(level=2, message="Configuration & Build Parameters:")
    utils.print_update(level=3, message=f"Year: {year}")
    utils.print_update(level=3, message=f"Copperplate Mode: {copperplate}")
    utils.print_update(level=3, message=f"Capacity Choice: {capacity_choice}")
    utils.print_update(level=3, message=f"TX Line Infinity: {tx_line_infinity}")
    utils.print_update(level=3, message=f"Include VRE investments: {include_vre_investments}")
    utils.print_update(level=3, message=f"Reservoir representation: {reservoir_representation}")
    utils.print_update(level=3, message=f"Water-use policy: {water_policy}")
    utils.print_update(level=3, message=f"Release multiplier: {release_multiplier}")
    utils.print_update(level=3, message=f"Solve network: {solve_network}")
    utils.print_update(level=3, message=f"Output File: {solved_network_save_to if solved_network_save_to else 'Default from config'}")

    valid_representations = {"A_full_cascade", "B_two_store", "C_single_bucket"}
    if reservoir_representation not in valid_representations:
        raise ValueError(
            f"Unsupported reservoir_representation={reservoir_representation!r}; "
            f"expected one of {sorted(valid_representations)}"
        )
    valid_water_policies = {
        "evidence_based",
        "static_minimum_only",
        "no_minimum_release",
        "release_uncertainty_lower_bound",
        "release_uncertainty_upper_bound",
    }
    if water_policy not in valid_water_policies:
        raise ValueError(
            f"Unsupported water_policy={water_policy!r}; expected one of "
            f"{sorted(valid_water_policies)}"
        )
    if not np.isfinite(release_multiplier) or release_multiplier < 0:
        raise ValueError("release_multiplier must be finite and non-negative")
    
    # Validate required data files
    utils.print_update(level=2, message="Validating required data files...")
    required_files = [
        pypsa_aparser.data_cfg["output"]["base_network"],
        pypsa_aparser.data_cfg["output"]["disaggregate_load"]["res_path"],
        pypsa_aparser.data_cfg["output"]["disaggregate_load"]["csmi_path"],
        pypsa_aparser.data_cfg["GADM"]["country_file_L1"] if copperplate else pypsa_aparser.data_cfg["GADM"]["country_file_L2"],
    ]
    
    missing = []
    for fpath in required_files:
        p = Path(fpath)
        if p.exists():
            size_kb = p.stat().st_size / 1024
            utils.print_update(level=3, message=f"✓ {fpath} ({size_kb:.1f} KB)")
        else:
            utils.print_update(level=3, message=f"✗ {fpath} (MISSING)")
            missing.append(fpath)
    
    if missing:
        raise FileNotFoundError(
            f"Missing {len(missing)} required data file(s):\n" + 
            "\n".join([f"  - {p}" for p in missing])
        )

    cfg=pypsa_aparser.data_cfg

    # (1) Load files
    utils.print_update(level=2,message='Creating PyPSA network...')
    # For PyPSA 0.28+, custom component attributes are handled differently
    # Create a simple network and extend it dynamically as needed
    network = pypsa.Network()
    
    utils.print_update(level=2,message='Loading the prepared data-file paths for PyPSA-BC')
    network_path = cfg['output']['base_network']
    hydro_ror_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["ror"]
    hydro_res_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["res"]
    hydro_ror_water_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["ror_water"]
    wind_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["wind"]
    pv_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["solar"]
    tpp_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["tpp"]
    ff_infra_path = cfg["output"]["pypsa_dict"]["folder"] + cfg["output"]["pypsa_dict"]["ff_infrastructure"]


    network.add("Carrier","AC", co2_emissions=0)
    network.add("Carrier","inflow", co2_emissions=0)
    network.add("Carrier","Water", co2_emissions=0)
    network.add("Carrier","DC", co2_emissions=0)
    network.add("Carrier","HDG", co2_emissions=0)
    
    network.import_from_csv_folder(network_path)
    utils.print_update(level=3,message=f'Loaded the network from : {network_path}')
    
    ror_dict = utils.read_pickle(hydro_ror_path)
    utils.print_update(level=3,message=f'Loaded the ROR data from : {hydro_ror_path}')
    
    res_dict = utils.read_pickle(hydro_res_path)
    utils.print_update(level=3,message=f'Loaded the Reservoir data from : {hydro_res_path}')
    
    ror_water_dict = utils.read_pickle(hydro_ror_water_path)
    utils.print_update(level=3,message=f'Loaded the Water inflow data from : {hydro_ror_water_path}')
    
    wind_dict = utils.read_pickle(wind_path)
    utils.print_update(level=3,message=f'Loaded the Wind resources data from : {wind_path}')
    
    pv_dict = utils.read_pickle(pv_path)
    utils.print_update(level=3,message=f'Loaded the Solar photovoltaic resources data from : {pv_path}')
    
    tpp_dict = utils.read_pickle(tpp_path)
    utils.print_update(level=3,message=f'Loaded the Thermal power resources data from : {tpp_path}')
    
    ff_infra_dict = utils.read_pickle(ff_infra_path)
    utils.print_update(level=3,message=f'Loaded the Fossil fuel resources data from : {ff_infra_path}')

    # (2) Set time-slicing
    utils.print_update(level=2,message='Updating ROR time-slices...')
    network.set_snapshots(ror_dict['BC_ABN_GSS']['p_max_pu'].index) # UPDATE REQUIRED

    # (3) Start adding assets
    utils.print_update(level=2,message='Adding power generation assets to the network...')
    add_hydro_ror_assets(network, ror_dict)
    add_hydro_res_assets(network, res_dict) # NOTE: None for AB
    add_hydro_ror_water_assets(network, ror_water_dict) # NOTE: None for AB
    if "minimum_release_m3_per_hour" in network.links.columns:
        static_minima = pd.to_numeric(
            network.links["minimum_release_m3_per_hour"], errors="coerce"
        ).fillna(0.0)
        if water_policy == "no_minimum_release":
            network.links.loc[:, "minimum_release_m3_per_hour"] = 0.0
        else:
            network.links.loc[:, "minimum_release_m3_per_hour"] = (
                static_minima * float(release_multiplier)
            )

    if release_schedule_path is None:
        release_schedule_path = Path(
            cfg.get("inventory", {}).get(
                "release_schedules",
                "studies/chapter2/inputs/release_schedules.csv",
            )
        )
    else:
        release_schedule_path = Path(release_schedule_path)
    schedule_policies = {
        "evidence_based",
        "release_uncertainty_lower_bound",
        "release_uncertainty_upper_bound",
    }
    if water_policy in schedule_policies and release_schedule_path.exists():
        release_schedules = pd.read_csv(release_schedule_path)
        release_schedules["minimum_release_m3_per_hour"] = (
            pd.to_numeric(
                release_schedules["minimum_release_m3_per_hour"], errors="raise"
            )
            * float(release_multiplier)
        )
        attached_routes = attach_route_release_schedules(network, release_schedules)
        if attached_routes:
            utils.print_update(
                level=2,
                message=(
                    f"Attached {attached_routes} evidence-labelled, route-specific "
                    f"release schedule(s) from {release_schedule_path}"
                ),
            )
    representation_summary = pd.DataFrame()
    if reservoir_representation != "A_full_cascade":
        topology_path = Path(
            cascade_topology_path
            or "studies/chapter2/inputs/cascade_evidence_register.csv"
        )
        if not topology_path.exists():
            raise FileNotFoundError(
                f"Cascade topology required by {reservoir_representation} is missing: "
                f"{topology_path}"
            )
        topology = pd.read_csv(topology_path, low_memory=False)
        representation_result = aggregate_cascade_representation(
            network,
            topology,
            reservoir_representation,
        )
        representation_summary = representation_result.summary
        utils.print_update(
            level=2,
            message=(
                f"Applied {reservoir_representation} across "
                f"{representation_summary['cascade_group'].nunique()} cascade groups; "
                "storage, inflows, turbine capacity, electrical buses, and external "
                "release routes were held fixed."
            ),
        )
    release_uncertainty_audit = pd.DataFrame()
    uncertainty_cases = {
        "release_uncertainty_lower_bound": "lower_bound",
        "release_uncertainty_upper_bound": "upper_bound",
    }
    if water_policy in uncertainty_cases:
        uncertainty_path = Path(
            release_uncertainty_protocol_path
            or "data/validation/chapter2/release_uncertainty_protocol.yaml"
        )
        if not uncertainty_path.exists():
            raise FileNotFoundError(
                f"Release uncertainty protocol is missing: {uncertainty_path}"
            )
        release_uncertainty_audit = apply_release_uncertainty_case(
            network,
            uncertainty_path,
            uncertainty_cases[water_policy],
        )
        utils.print_update(
            level=2,
            message=(
                f"Applied scenario-labelled release uncertainty policy to "
                f"{len(release_uncertainty_audit)} unresolved stations; these "
                "bounds are not observed or regulatory release rules."
            ),
        )
    add_pv_assets(network, pv_dict)
    add_wind_assets(network, wind_dict)
    add_ff_infra(network, ff_infra_dict) # NOTE: AB Additon 
    add_tpp_assets(network, tpp_dict)
    # add_cogen_assets(network, cogen_dict)
    
    # (4) add carriers outside of default ELC
    
    # network.add("Carrier","NG", co2_emissions=1.0)
    # network.add(class_name = 'Store',
    #             name ="Global NG Store",
    #             bus = "Global NG Bus",
    #             e_nom = 1e10,
    #             e_initial=1e10)
    
    # (5) Clustering
    
    # (6)  get busmap dictionary
    # NOTE: Centroid calculation will need an update and validation.
    # bus_dict = {name:0 for name in pd.read_csv(bus_path)['name'].tolist()}
    # NOTE: GADM switching to BC entire province here just for time being
    
    # region =  cfg["output"]["build_model"]["region_res"] # "single": BC as a single region, "multiple": BC split into 28 regional districts 
    
    if not copperplate:
        utils.print_update(level=2,message="Preparing network for multi regional nodes...")
        # region == 'multiple':
        # # geo_file = cfg["data"]["gadm"]["bc"] #"/mnt/c/Users/pmcw9/Delta-E/PICS/Data/regions/gadm41_CAN_2.json"
        # geo_file = r'/home/pmcwhannel/repos/PyPSA_BC/data/regions/AESO-Planning-Areas-2020-06-23'
        # regional_gdf = gpd.read_file(geo_file)

        # # (7A) Add new gadm regions as buses

        # Fetch or load GADM level 2 data for Canada
        regional_boundaries = Path(cfg["GADM"]["country_file_L2"])
        regional_boundaries.parent.mkdir(parents=True, exist_ok=True)
        
        if regional_boundaries.exists():
            gdf = gpd.read_file(regional_boundaries)
        else:
            import pygadm
            gdf = pygadm.Items(name="Canada", content_level=2)
            if gdf.crs is None:
                gdf = gdf.set_crs("EPSG:4326")
            gdf.to_file(regional_boundaries, driver="GeoJSON")
        
        # Debug: Check available columns and BC entries
        utils.print_update(level=3, message=f"GADM columns: {list(gdf.columns)}")
        utils.print_update(level=3, message=f"Unique NAME_1 values: {gdf['NAME_1'].unique()[:5]}")  # Show first 5
        
        gadm_bc = gdf[gdf["NAME_1"] == "British Columbia"]
        if len(gadm_bc) == 0:
            # Try alternative names
            gadm_bc = gdf[gdf["NAME_1"].str.contains("British", case=False, na=False)]
            if len(gadm_bc) == 0:
                utils.print_update(level=1, message="Warning: No BC found with 'British Columbia' or 'British' pattern")
                gadm_bc = gdf  # Use all data
        busmap_dict = get_bc_busmap_dict(network, gadm_bc)

        # (7A) Add new gadm regions as buses
        utils.print_update(level=3,message="Creating regional buses...")
        create_adm_buses(network, gadm_bc, busmap_dict)

    else:
        utils.print_update(level=3,message="Preparing network for single node (aggregated regions)")
        # Fetch or load GADM level 1 data for Canada (provincial boundaries)
        utils.print_update(level=3, message="Loading GADM level 1 administrative boundaries...")
        geojson_file = cfg["GADM"]["country_file_L1"]
        geojson_file_path = Path(geojson_file)
        geojson_file_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            if geojson_file_path.exists():
                gdf = gpd.read_file(geojson_file)
                utils.print_update(level=3, message=f"✓ Loaded GADM L1 from file: {geojson_file}")
            else:
                # Fetch from pygadm and save
                import pygadm
                utils.print_update(level=3, message="Fetching GADM L1 data from pygadm...")
                gdf = pygadm.Items(name="Canada", content_level=1)
                if gdf.crs is None:
                    gdf = gdf.set_crs("EPSG:4326")
                gdf.to_file(geojson_file, driver="GeoJSON")
                utils.print_update(level=3, message=f"✓ Saved GADM L1 to: {geojson_file}")
            
            mask = gdf["NAME_1"] == "BritishColumbia"
            busmap_dict = get_single_region_busmap_dict(network, gdf.loc[mask,:])
        except Exception as e:
            utils.print_update(level=1, message=f"Error loading GADM L1 data: {e}")
            raise

        # (7A) Add new gadm regions as buses: (In this case just the single BC bus)
        network.add(class_name="Bus",
                    name="BC",
                    x=gdf.loc[mask,:].geometry.centroid.x.iloc[0],
                    y=gdf.loc[mask,:].geometry.centroid.y.iloc[0],
                    v_nom = 300 # voltage assumed
                    )
   

    # (7B) Add new trade buses
    # NOTE: trade buses have to be added to busmap_dict after the administrative buses are added atm.
    # busmap_dict["500_BCUS1_INT"] = "US Trade"
    # busmap_dict["230_BCUS2_INT"] = "US Trade"
    # busmap_dict["138_BCAB1_IPT"] = "AB Trade"
    # busmap_dict["500_BCAB2_IPT"] = "AB Trade"
    # busmap_dict["138_BCAB3_IPT"] = "AB Trade"
    # busmap_dict["138_BCAB4_IPT"] = "AB Trade"

    # create_trade_buses(network)


    # NOTE: Think about how to aggregate the new lines...
    # network.lines['s_nom'] = 50000 # 3700 no good. Good at 3800. Good at 4000.

    # (8) Add load (BC)

    (start_time,end_time)=pypsa_aparser.snapshot
    utils.print_update(level=2,message= f"Snapshot loaded: {start_time}-{end_time}")

    utils.print_update(level=2,message= "Checking loads data")
    res_load = pd.read_csv(cfg['output']['disaggregate_load']['res_path'],
                            index_col=0, parse_dates=True).loc[start_time:end_time]
    utils.print_update(level=3,message= f"Residential load data loaded from: {cfg['output']['disaggregate_load']['res_path']} ")

    csmi_load = pd.read_csv(cfg['output']['disaggregate_load']['csmi_path'],
                            index_col=0, parse_dates=True).loc[start_time:end_time]
    utils.print_update(level=3,message= f"Commercial, Small and Medium Industries load data loaded from: {cfg['output']['disaggregate_load']['res_path']} ")
    
    # # NOTE: Update coming for Stikine and CentralCoast.z
    # NOTE: May be issues here with assignment of new assets
    # These regions have no way to serve load or do not have load according to CEEI
    utils.print_update(level=2,message="Calibrating load for network buses...")
    if not copperplate:
        for col in res_load.columns: 
            if col not in ["Stikine", "CentralCoast", "NorthernRockies"]: 
                load_ts = res_load[col] + csmi_load[col]
                network.add("Load", "{} ELC Load".format(col), bus=col, p_set=load_ts) # Make load much smaller
                # network.add(class_name="Generator",
                #     name='Backstop {}'.format(col),
                #     bus=col,
                #     p_nom=10000.0,
                #     marginal_cost=1000,
                #     )
            else:
                pass
    else:
        total_load_ts = res_load.sum(axis=1) + csmi_load.sum(axis=1)
        network.add("Load", "BC ELC Load", bus="BC", p_set=total_load_ts) # Make load much smaller



    # (9) Determine all components to relink
    utils.print_update(level=2,message="Refactoring,aggregating and mapping buses to regions...")
    replace_bus_refs(network.generators,'bus', busmap_dict)
    replace_bus_refs(network.links,'bus0', busmap_dict)
    replace_bus_refs(network.links,'bus1', busmap_dict)
    replace_bus_refs(network.links,'bus2', busmap_dict)
    replace_bus_refs(network.lines,'bus0', busmap_dict)
    replace_bus_refs(network.lines,'bus1', busmap_dict)
    # replace_bus_refs(network.loads,'bus', busmap_dict)

    # (10) Remove the old components which are no longer needed (lines, buses)
    utils.print_update(level=2,message="Cleaning the non-aggregated lines and buses that are not required for solving the network...")
    remove_old_components(network, busmap_dict)

    #################### ADDING ALL NEW ASSETS BEGINS HERE ##############
    # (11) Add VRE Expansion
    # capex: $M-CAD / (potential capacity)
    # potential_capacity: MW
    # CF_mean: capacity factor
    # p_lcoe: MW-hr / $M-CAD-per-MW-Installed
    # load pv/solar, wind, and battery assets?
    resource_options_data = Path('data/processed_data')

    if include_vre_investments:
        utils.print_update(level=2,message="Loading Resource Options (non-existing future resources)...")
        solar_resources=resource_options_data/'solar/potential'
        wind_resources=resource_options_data/'wind/potential'

        pv_sites = pd.read_csv(solar_resources/f'resource_options_investments_solar_{year}.csv',index_col='cluster_id')# utils.read_pickle(os.path.join(temp_folder,cfg_complete['results']['linking']['clusters_topSites']['solar']))
        pv_ts = pd.read_csv(solar_resources/'resource_options_solar_timeseries.csv',index_col='time', parse_dates=True)# utils.read_pickle(os.path.join(temp_folder,cfg_complete['results']['linking']['clusters_CFts_topSites']['solar']))
        utils.print_update(level=3,message=f"Solar sites and profiles loaded from {solar_resources} ")
   
        wind_sites = pd.read_csv(wind_resources/f'resource_options_investments_wind_{year}.csv',index_col='cluster_id')# utils.read_pickle(os.path.join(temp_folder,cfg_complete['results']['linking']['clusters_topSites']['wind']))
        wind_ts = pd.read_csv(wind_resources/'resource_options_wind_timeseries.csv',index_col='time', parse_dates=True)# utils.read_pickle(os.path.join(temp_folder,cfg_complete['results']['linking']['clusters_CFts_topSites']['wind']))
        utils.print_update(level=3,message=f"Wind sites and profiles loaded from {wind_resources} ")
    

        add_vre_expansion_sites(network, pv_sites, pv_ts, vre_type='Solar',capacity_choice=capacity_choice)
        add_vre_expansion_sites(network, wind_sites, wind_ts, vre_type='Wind',capacity_choice=capacity_choice)
    else:
        utils.print_update(level=2, message="Skipping VRE investment sites for this build run.")
    
    committed_sites_solar=pd.read_csv('data/processed_data/solar/committed/BCH_CFP24_solar.csv',index_col='project_name')
    committed_sites_solar_ts=pd.read_csv('data/processed_data/solar/committed/BCH_CFP24_solar_ts.csv',index_col='time',parse_dates=True)
    committed_sites_wind=pd.read_csv('data/processed_data/wind/committed/BCH_CFP24_wind.csv',index_col='project_name')
    committed_sites_wind_ts=pd.read_csv('data/processed_data/wind/committed/BCH_CFP24_wind_ts.csv',index_col='time',parse_dates=True)
    
    committed_sites_solar_year=committed_sites_solar.iloc[committed_sites_solar['start_year'].values <= year]
    committed_sites_wind_year=committed_sites_wind.iloc[committed_sites_wind['start_year'].values <= year]
    
    if committed_sites_solar_year.empty:
        utils.print_update(level=2,message=f"No Committed Solar sites for {year}")
    else:
        add_vre_committed_sites(network, committed_sites_solar, committed_sites_solar_ts, vre_type='Solar')
        
    if committed_sites_wind_year.empty:
        utils.print_update(level=2,message=f"No Committed Wind sites for {year}")
    else:
        add_vre_committed_sites(network, committed_sites_wind, committed_sites_wind_ts, vre_type='Wind')
        
    
    utils.print_update(level=3,message="Resources options added to pypsa network.")
    
    ### Add Nuclear Resource Options
    if include_vre_investments:
        utils.print_update(level=3,message="Loading Nuclear Resources")
        nuc_resources=resource_options_data/'nuclear/potential'
        
        nuc_sites = pd.read_csv(nuc_resources/f'resource_options_investments_nuc_{year}.csv',index_col='name')# utils.read_pickle(os.path.join(temp_folder,cfg_complete['results']['linking']['clusters_topSites']['solar']))
        utils.print_update(level=3,message=f"Nuclear site(s) loaded from {nuc_resources} ")

        add_nuc_expansion_sites(network, nuc_sites, capacity_choice=capacity_choice, bus_fix='GreaterVancouver')
    else:
        utils.print_update(level=2, message="Skipping nuclear investment sites for this build run.")

# data/processed_data/nuclear/potential/resource_options_investments_nuc_2032.csv
    
    # # (12) Add trade load
    # # NOTE: Needs to be updated in line with what is done in BC_Nexus
    # # add_trade(network, cfg)

    # NOTE: Temporary fix on ramping
    # network.generators["ramp_limit_down"] = 1.0
    # network.generators["ramp_limit_up"] = 1.0
    # network.links["ramp_limit_down"] = 1.0
    # network.links["ramp_limit_up"] = 1.0
    
    # Add Backstops
    # NOTE: Watch out for when cogen is added here....
    utils.print_update(level=100,message="Creating BACKSTOP generators to the network...")
    for bus in network.loads.bus.unique():
        network.add(class_name="Generator",
                    name='Backstop {}'.format(bus),
                    type='backstop',
                    bus=bus,
                    p_nom=50000,
                    marginal_cost=1000,
                    p_nom_extendable=False,
                    capital_cost=999999
                    )
                
    # Aggregate lines
    utils.print_update(level=2,message="Aggregating inter-zones lines...")
    aggregate_lines(network)
    repair_bus_references(network)
    network.sanitize()
    coerce_network_string_labels(network)
    prepare_and_validate_line_parameters(network)
    
    if tx_line_infinity:
        network=make_unlimited_line_capacity(network)
    
    # Model construction and optimization are separate workflow stages.  Keep
    # the legacy default (solve_network=True) for direct callers, while the
    # Snakemake base-model workflow exports an unsolved network first.
    if solve_network:
        utils.print_update(level=2,message="Optimizing the built network...")
        solver_name = _choose_solver_name()
        utils.print_update(level=3, message=f"Using solver: {solver_name}")
        status, termination_condition = network.optimize(
            solver_name=solver_name,
            extra_functionality=add_cascade_constraints,
        )
        utils.print_update(
            level=2,
            message=(
                f"Optimization completed with status: {status}, "
                f"termination: {termination_condition}"
            ),
        )
    else:
        status, termination_condition = "NOT_RUN", "model_preparation_only"
        utils.print_update(
            level=2,
            message="Model preparation complete; optimization intentionally skipped.",
        )
    # Note: In PyPSA 0.28+, network.model is automatically created during optimize()
    # No need to call create_model() separately

    # Save network
    if solved_network_save_to is None:
        base_output = Path(cfg["output"]["build_model"]["fname"])
        if base_output.suffix.lower() == ".nc":
            solved_network_save_to = base_output.with_name(f"{base_output.stem}_{year}.nc")
        else:
            solved_network_save_to = base_output.with_name(f"{base_output.name}_{year}.nc")
    
    solved_network_save_to.parent.mkdir(exist_ok=True,parents=True)
    network.export_to_netcdf(solved_network_save_to)
    network_state = "Solved" if solve_network else "Prepared (unsolved)"
    utils.print_update(level=1,message=f"{network_state} network saved to: {solved_network_save_to}")
    
    # Generate markdown report
    utils.print_update(level=2,message="Generating build model report...")
    report_path = (
        Path(build_report_save_to)
        if build_report_save_to is not None
        else Path("reports") / f"build_model_report_{year}.md"
    )
    report_path.parent.mkdir(exist_ok=True, parents=True)
    if not representation_summary.empty:
        representation_summary.to_csv(
            report_path.with_name(f"{report_path.stem}_hydraulic_aggregation.csv"),
            index=False,
        )
    if not release_uncertainty_audit.empty:
        release_uncertainty_audit.to_csv(
            report_path.with_name(f"{report_path.stem}_release_uncertainty.csv"),
            index=False,
        )
    
    report_content = f"""# PyPSA-BC Build Model Report
    
**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Model Configuration

| Parameter | Value |
|-----------|-------|
| Year | {year} |
| Copperplate Mode | {copperplate} |
| Capacity Choice | {capacity_choice} |
| TX Line Infinity | {tx_line_infinity} |
| Reservoir Representation | {reservoir_representation} |
| Water-use Policy | {water_policy} |
| Release Multiplier | {release_multiplier} |
| Optimization requested | {solve_network} |
| Release Evidence Boundary | {'Scenario sensitivity; not observed/legal' if not release_uncertainty_audit.empty else 'No uncertainty envelope applied'} |
| Output File | {solved_network_save_to} |

## Network Summary

| Component | Count |
|-----------|-------|
| Buses | {len(network.buses)} |
| Lines | {len(network.lines)} |
| Transformers | {len(network.transformers)} |
| Generators | {len(network.generators)} |
| Stores | {len(network.stores)} |
| Links | {len(network.links)} |
| Time Snapshots | {len(network.snapshots)} |

## Data Sources

### GADM Geographic Data
- **Level**: {'2 (Regional Districts)' if not copperplate else '1 (Provinces)'}
- **Source**: pygadm or cached/local GADM files
- **Cached Locations**: 
  - Level 1: `data/downloaded_data/GADM/gadm_can_l1.geojson`
  - Level 2: `data/downloaded_data/GADM/gadm_can_l2.geojson`

### PyPSA Assets
- **Base Network**: {cfg['output']['base_network']}
- **Hydro RoR**: {cfg['output']['pypsa_dict']['folder']}{cfg['output']['pypsa_dict']['ror']}
- **Hydro Reservoirs**: {cfg['output']['pypsa_dict']['folder']}{cfg['output']['pypsa_dict']['res']}
- **Wind**: {cfg['output']['pypsa_dict']['folder']}{cfg['output']['pypsa_dict']['wind']}
- **Solar**: {cfg['output']['pypsa_dict']['folder']}{cfg['output']['pypsa_dict']['solar']}
- **Thermal**: {cfg['output']['pypsa_dict']['folder']}{cfg['output']['pypsa_dict']['tpp']}

## Optimization Status

- **Status**: {status}
- **Termination Condition**: {termination_condition}

---
*Report generated by PyPSA-BC workflow*
"""
    
    with open(report_path, 'w') as f:
        f.write(report_content)
    
    utils.print_update(level=1, message=f"Build model report saved to: {report_path}")

if __name__ == '__main__':
    main()
