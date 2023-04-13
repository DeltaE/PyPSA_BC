import pypsa
from bc_power import utils, hydro
import json
import pandas as pd


def is_terminal_stage(down_rid):
    '''
    This function checks whether a cascade is a terminal cascade or not
    return (boolean): True if stage is terminal otherwise False.
    '''
    terminal_list = {"DEFAULT","Columbia (U.S.)", "Brilliant", "Lower Bonnington",
                   "Upper Bonnington"}
    if down_rid in terminal_list:
        return True
    elif down_rid:
        return 
    else:
        return False


def get_reservoir_dict(site, reservoir, inflow, res_list):
    '''
    This function is used to return a dictionary with 7 components to be used to build
    a hydroelectric reservoirs in PyPSA.
    Components:
    1) Water bus: Bus connected to discharge, spill, and connected by links to reservoir bus. (No storage on water bus)
    2) Reservoir bus: Bus directly connected to reservoir storage.
    3) Store link: Link used to move water from water bus into the reservoirs storage
    4) Release link: Link used to release water from the reservoir bus and thereby storage into the water bus for use.
    5) Reservoir store: Store used to store water for the reservoir being modelled.
    6) Inflow generator: Generator used to be source of inflow into the system.
    7) Discharge link: Link used to produces power and then release the same water to next stage.
    8) Spill link: Link use to spill water to the next stage.
    
    '''
    # IDs
    aid = site.name
    up_rid = site["upper_reservoir_id"]
    down_rid = site["lower_reservoir_id"]

    # parameters
    # NOTE: Unit change should happen here if desired otherwise leave as is..
    max_storage = reservoir["max_storage"] - reservoir["min_storage"]
    max_inflow = max(inflow)


    res_dict = {}
    # 1) add water_bus
    # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
    carrier = "Water"
    class_name = "Bus"
    res_dict['water bus'] = {"class_name":class_name,
                                "name":" ".join([aid,carrier,class_name]),
                                "carrier":carrier,
                                }

    # 2) add reservoir_bus
    # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
    carrier = "Water"
    class_name = "Bus"
    res_dict['reservoir bus'] = {"class_name":class_name,
                                    "name":" ".join([up_rid,carrier,class_name]),
                                    "carrier":carrier,
                                    }

    # 3) add store_link
    # name: {ASSET_ID} Store {CLASS_NAME}
    class_name = "Link"
    res_dict['store link'] = {"class_name":"Link",
                                    "name": " ".join([aid,"Store",class_name]),
                                    "bus0": res_dict['water bus']['name'],
                                    "bus1": res_dict['reservoir bus']['name'],
                                    "efficiency":1., # mass balance
                                    "p_nom":max_storage, # Should be adjusted based on timestep eventually
                                    } 
    
    # 4) get release_link
    # name: {ASSET_ID} Release {CLASS_NAME}
    class_name = "Link"
    res_dict['release link'] = {"class_name":class_name,
                                    "name": " ".join([aid,"Release",class_name]),
                                    "bus0": res_dict['reservoir bus']['name'],
                                    "bus1": res_dict['water bus']['name'],
                                    "efficiency":1., # mass_balance
                                    "p_nom":max_storage + max_inflow,
                                    } 
    

    # 5) get reservoir store
    # name: {UPPER_RESERVOIR_ID} Store {CLASS_NAME}
    class_name = "Store"
    res_dict['reservoir store'] = {"class_name":class_name,
                                        "name":" ".join([up_rid,"Reservoir",class_name]),
                                        "bus":res_dict['reservoir bus']['name'],
                                        "e_nom":max_storage
                                        }
    
    # 6) add inflow generator
    # name: {UPPER_RESERVOIR_ID} Inflow {CLASS_NAME}
    class_name = "Generator"
    res_dict['inflow generator'] = {"class_name":class_name,
                                        "name": " ".join([up_rid,"Inflow",class_name]),
                                        "bus": res_dict['reservoir bus']['name'],
                                        "carrier": "inflow",
                                        "efficiency":1., # mass_balance
                                        "p_nom":max_inflow, # max(inflow series)
                                        "p_set":inflow,
                                        "p_max_pu":[flow / max_inflow if flow != 0 else 0 for flow in inflow],
                                        "p_min_pu":[flow / max_inflow if flow != 0 else 0 for flow in inflow],
                                        }
    

    # Terminal check
    class_name = "Bus"
    carrier = "Water"
    if down_rid in res_list: # true if not terminal
        # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
        downstream_bus = " ".join([down_rid, carrier, class_name])
    else:
        # name: {ASSET_ID} {CARRIER} {CLASS_NAME}
        cascade_name = site["cascade_group"]
        downstream_bus = " ".join([cascade_name,carrier,class_name])
        # Add bus for the terminal reservoir
        res_dict['terminal bus'] = {"class_name":"Bus",
                                    "name":downstream_bus,
                                    "carrier":"water",
                                    }
        
        # Add spill store
        res_dict['terminal store'] = {"class_name":"Store",
                                        "name":cascade_name,
                                        "bus":downstream_bus,
                                        "e_nom":1e45 # The max storage needs to retain all possible water in the model horizon
                                        }

    # 7) get discharge link
    elc_bus_name = " ".join([site["connecting_node_code"],"ELC","Bus"]) # ELC BUS NAME based on existing schema
    q_rated = float(site['max_water_discharge']) # This flow would have to be modified for non-hourly models
    res_dict['discharge link'] = {"class_name":"Link",
                                        "name": " ".join([aid,"discharge link"]),
                                        "bus0": res_dict['water bus']['name'],
                                        "bus1": elc_bus_name, # NEED TO CREATE connection to ELC (ELC bus name?)
                                        "bus2": downstream_bus,
                                        "marginal_cost":site["variable_om_cost_USD_per_MWh"],
                                        "efficiency":site['capacity'] / q_rated,
                                        "efficiency2":1., # mass balance
                                        "p_nom":q_rated, # Should be derived to ensure larger than max(inflow, spill + discharge)
                                        } 
    # 8) get spill link
    if site['max_spill'] == "DEFAULT":
        spill = 9999 # Fill a value
    else:
        spill = float(site['max_spill'])

    if spill == 0: # Later this should be replaced to allow 0 spill
        spill = 9999 # Fill a value

    res_dict['spill_link'] = {"class_name":"Link",
                                "name": " ".join([aid,"spill link"]),
                                "bus0": res_dict['water bus']['name'],
                                "bus1": downstream_bus,
                                "efficiency":1., # mass_balance
                                "p_nom":spill, # Should be derived to ensure larger than max(inflow, spill + discharge)
                                }
    
    return res_dict

def get_ror_dict(site, ror_ts):
    '''
    This function creates the dictionary for a single RoR facility.
    '''
    # name: {ASSET_ID} RoR Generator
    # bus" {CONNECTING_NODE_CODE} ELC Bus
    name = " ".join([site["asset_id"], "RoR", "Generator"])
    elc_bus = " ".join([site["connecting_node_code"],"ELC", "Bus"])

    return {"class_name":"Generator",
            "name":name,
            "bus":elc_bus,
            "p_nom":site['capacity'],
            "marginal_cost":site["variable_om_cost_USD_per_MWh"],
            "p_nom_extendable":False, # Site already built
            # "capital_cost":site[], # no applicable since built
            "p_max_pu":ror_ts}



def write_reservoir_dict(hydro_sites, hydro_res, res_inflows, cfg):
    '''
    This function creates and writes the reservoir dictionary which contains all the necessary components
    to be read in by pypsa to instantiate the reservoirs in pypsa.
    '''
    subset = ["asset_id","latitude","longitude",]
    sum_list = ["capacity", "annual_avg_energy", "ramp_up", "ramp_down"] # list of parameters to aggregate (assume ramping applies per turbine and asset)
    temp_df = hydro_sites[hydro_sites["hydro_type"].str.contains("reservoir")] # filer for reservour hydro assets
    agg_hydro_sites = temp_df.groupby(by="asset_id", group_keys=False).apply(lambda x:
                                                    hydro.merge_assets(x, subset, sum_list))
    # Make sure to aggregate the assets

    # Create dictionary for json
    # Asset_ID -> components -> attributes
    res_dict = {}
    res_list = hydro_res['asset_id'].tolist() # list of unique reservoirs modelled here
    for _,site in agg_hydro_sites.iterrows():
        # IDs needed to index inflow and reservoir
        aid = site.name
        up_rid = site["upper_reservoir_id"]

        # find matching reservoirs
        reservoir = hydro_res[hydro_res["asset_id"] == up_rid].iloc[0]
        inflow = res_inflows[up_rid]
        res_dict[aid] = get_reservoir_dict(site, reservoir, inflow, res_list) 

    utils.write_pickle(res_dict, cfg["pypsa_dict"]["hydro_res"])

def write_ror_dict(hydro_sites, ror_series, cfg):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing RoR facilities in PyPSA.
    '''
    ror_dict = {}
    temp_df = hydro_sites[hydro_sites["hydro_type"].str.contains("ror")]
    for _,site in temp_df.iterrows():
        aid = site["asset_id"]
        ror_ts = ror_series[aid]
        ror_dict[aid] = get_ror_dict(site, ror_ts)


    # write json now
    utils.write_pickle(ror_dict, cfg["pypsa_dict"]["hydro_ror"])

def main():
    '''
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
    '''
    # Read in configuration file
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    hydro_sites = pd.read_csv(cfg["hydro_prep"]["hydro_generation"])
    hydro_res = pd.read_csv(cfg["hydro_prep"]["hydro_reservoir"]) # Purely reservoir information
    res_inflows = pd.read_csv(cfg["reservoir_inflows"]["output"])
    ror_series = pd.read_csv(cfg["ror_inflows"]["ror_outfile"])

    # (1) Write json dictionaries for the reservoirs.
    write_reservoir_dict(hydro_sites, hydro_res, res_inflows, cfg)

    # (2) Write json dictionaries for the RoR facilities.
    write_ror_dict(hydro_sites, ror_series, cfg)
    # save res_dict to a json
    
if __name__ == '__main__':
    main()