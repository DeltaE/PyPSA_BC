import pypsa
from bc_power import utils, hydro
import json
import pandas as pd


def get_vre_params(gen_generic, cfg):
    '''
    This function gets generic wind/solar parameters. Currently, it only pull form CODERS
    and return the onshore wind information (currently only costs are used). Will be updated in the future.
    '''

    if cfg['vre']['type'] == 'wind':
        return gen_generic[gen_generic["generation_type"] == "Wind_onshore"]
    elif cfg['vre']['type'] == 'solar':
        return gen_generic[gen_generic["generation_type"] == "Solar_PV"]
    else:
        vre_type = cfg['vre']['type'] 
        print(f"error: {vre_type} is not implemented yet!")
        exit(1) 



def get_vre_dict(site, site_ts, vre_params, bus_dict, cfg): # site, site_ts, vre_params, bus_dict, cfg
    '''
    This function cre
    '''
    # name: {ASSET_ID} Wind Generator
    # bus" {CONNECTING_NODE_CODE} ELC Bus
    name = " ".join([site["asset_id"], cfg['vre']['type'].title(), "Generator"])
    elc_bus = utils.get_gen_bus(site["connecting_node_code"], bus_dict)
    
    p_nom = site['Install capacity']
    
    return {"class_name":"Generator",
            "name":name,
            "bus":elc_bus,
            "p_nom":p_nom,
            "marginal_cost":0.000001, # vre_params["variable_om_cost_USD_per_MWh"] -> not right type
            "p_nom_extendable":False, # Site already built
            # "capital_cost":site[], # no applicable since built
            "p_max_pu":site_ts.apply(lambda x: min(x/ site['Install capacity'],1))} # Needs to be renormalized to p_nom

def write_vre_dict(vre_assets, vre_ts, vre_params, bus_dict, cfg):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing vre facilities in PyPSA.
    '''
    vre_dict = {}
    # reduce to single assets
    vre_unique = vre_assets.drop_duplicates(subset=['asset_id', 'latitude','longitude'],keep='first')

    for _,site in vre_unique.iterrows():
        aid = site["asset_id"]
        site_ts = vre_ts[aid] # index timeseries
        vre_dict[aid] = get_vre_dict(site, site_ts, vre_params, bus_dict, cfg)

    # write pickle
    out_file = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"][cfg['vre']['type']]
    utils.write_pickle(vre_dict, out_file)

def main():
    '''
    This script takes creates the VRE dictionary for wind or solar assets for PyPSA_BC
    '''
    # Read in configuration file
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    start_time = cfg['scope']['temporal']['start'] 
    end_time = cfg['scope']['temporal']['end']

    gen_generic = pd.read_csv(cfg["coders"]["gen_generic"])
    vre_assets = pd.read_csv(cfg["vre"]['asset_path'])
    vre_ts = pd.read_csv(cfg['vre']['ts_path'], index_col=0, parse_dates=True).loc[start_time:end_time]
    buses = pd.read_csv(cfg["network"]["folder"] + "/buses.csv")['name'].tolist()

    # (0A) Create folders if they have not been created already
    utils.create_folder(cfg["pypsa_dict"]['components'])

    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    bus_dict = utils.create_standard_gen_bus_map(buses)


    # (0C) Get vre cost information
    vre_params = get_vre_params(gen_generic, cfg)
                                
    # (1) Write pickle dictionaries for the vre assets.
    write_vre_dict(vre_assets, vre_ts, vre_params, bus_dict, cfg)

    
if __name__ == '__main__':
    main()