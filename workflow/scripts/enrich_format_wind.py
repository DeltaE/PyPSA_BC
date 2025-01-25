import pypsa
from bc_combined_modelling import utils, hydro
import json
import pandas as pd
import sys

def get_vre_params(gen_generic, vre_selection):
    '''
    This function gets generic wind/solar parameters. Currently, it only pull fromm CODERS
    and return the onshore wind information (currently only costs are used). Can be updated in the future.
    '''

    if vre_selection == 'wind':
        return gen_generic[gen_generic["generation_type"] == "Wind_onshore"]
    elif vre_selection == 'solar':
        return gen_generic[gen_generic["generation_type"] == "Solar_PV"]
    else:
        vre_type = vre_selection
        print(f"error: {vre_type} is not implemented yet!")
        exit(1) 



def get_vre_dict(site, site_ts, bus_dict, vre_selection): # site, site_ts, vre_params, bus_dict, cfg
    '''
    This function takes in a VRE site (wind or pv) and creates a dictionary which can be used in PyPSA
    to add a generator which represents the corresponding VRE asset.
    site: Row of a Dataframe containing most of the parameters relevant to the asset of interest.
    site_ts: Timeseries of the generation of the wind asset 
    vre_params: This can be passed and would contain financial information for the asset.
    '''
    # name: {ASSET_ID} Wind Generator
    # bus" {CONNECTING_NODE_CODE} ELC Bus
    name = " ".join([site["asset_id"], vre_selection.title(), "Generator"])
    elc_bus = utils.get_gen_bus(site["connecting_node_code"], bus_dict)
    
    p_nom = site['Install capacity']
    
     # NOTE: Updated "marginal_cost" from hard coded to a value from gen generic
     # for onsshore wind (i.e. 0.0001 -> 0.00). Sometimes 0 marginal cost can cause error in Optimization.

    return {"class_name":"Generator",
            "name":name,
            "bus":elc_bus,
            "p_nom":p_nom,
            "marginal_cost":site["variable_om_cost_CAD_per_MWh"],
            "p_nom_extendable":False, # Site already built
            # "capital_cost":site[], # not applicable since built
            "p_max_pu":site_ts.apply(lambda x: min(x / p_nom,1))} # Needs to be renormalized to p_nom

def write_vre_dict(vre_assets, vre_ts, bus_dict, vre_selection, vre_path):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing vre facilities in PyPSA.
    '''
    vre_dict = {}
    # reduce to single assets
    # NOTE: 2023-10-05 Updated this to aggregate common asset_id rather than retain first
    # Retaining the first was not the correct approach since capacity was lost and unaccounted for.
    # vre_unique = vre_assets.drop_duplicates(subset=['asset_id', 'latitude','longitude'], keep='first')
    vre_unique = vre_assets.groupby(by="asset_id").agg({
        'Install capacity':'sum', **{col: 'first' for col in vre_assets.columns if col != 'Install capacity'}
        }).reset_index(drop=True)

    for _,site in vre_unique.iterrows():
        aid = site["asset_id"]
        site_ts = vre_ts[aid] # index timeseries
        vre_dict[aid] = get_vre_dict(site, site_ts, bus_dict, vre_selection)

    # write pickle
    out_file = vre_path
    utils.write_pickle(vre_dict, out_file)

def main():
    '''
    This script takes creates the VRE dictionary for wind or solar assets for PyPSA_BC
    '''
    # Read in configuration file
    config_file = r"config/config.yaml"   
    cfg_complete = utils.load_config(config_file)
    cfg=cfg_complete['pypsa']
    print(">>> formatting wind dataset initiates...")

    # start_time = cfg['params']['start']
    # end_time = cfg['params']['end']
    
    start_time = cfg_complete['province_mapping']['BC']['snapshots_tz_BC']['start'][0]
    end_time = cfg_complete['province_mapping']['BC']['snapshots_tz_BC']['end'][0]

    # These all derived from configuration file
    asset_path = cfg['output']['create_ext_wind_assets']['fname']
    ts_path = cfg['output']['create_ext_wind_ts']['fname']
    vre_selection = cfg['output']['enrich_format_wind']['vre_sel']
    vre_path = "data/pypsa_data/" + cfg['output']['pypsa_dict']['wind']


    # gen_generic = pd.read_csv(cfg['data']["coders"]["gen_generic"])
    
    vre_assets = pd.read_csv(asset_path)
    vre_ts = pd.read_csv(ts_path, index_col=0, parse_dates=True).loc[start_time:end_time]
    buses = pd.read_csv(cfg['output']['prepare_base_network']['folder'] + "/buses.csv")['name'].tolist()

    # (0A) Create folders if they have not been created already
    utils.create_folder(cfg['output']["pypsa_dict"]['folder'])

    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    # buses.append('138_604S_GSS') #NOTE: Need to look into this at somepoint.. (Why missing in buses.)
    bus_dict = utils.create_standard_gen_bus_map(buses)


    # # (0C) Get vre cost information
    # vre_params = get_vre_params(gen_generic, vre_selection)
                                
    # (1) Write pickle dictionaries for the vre assets.
    province = cfg['output']['prepare_base_network']['regions'][0]
    write_vre_dict(vre_assets, vre_ts, bus_dict, vre_selection, vre_path)
    print(">>> formatting wind dataset completed !")
    
if __name__ == '__main__':
    # 
    main()
    # main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    # "data/pypsa_data/bc_solar_assets.csv", "results/interim/bc_solar_ts.csv", "solar", "/home/pmcwhannel/repos/PyPSA_BC/results/pypsa-components/bc_solar.pickle"