import pypsa
from bc_power import utils, hydro
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



def get_vre_dict(site, site_ts, vre_params, bus_dict, vre_selection): # site, site_ts, vre_params, bus_dict, cfg
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
            "marginal_cost":vre_params["variable_om_cost_USD_per_MWh"].iloc[0],
            "p_nom_extendable":False, # Site already built
            # "capital_cost":site[], # not applicable since built
            "p_max_pu":site_ts.apply(lambda x: min(x / p_nom,1))} # Needs to be renormalized to p_nom

def write_vre_dict(vre_assets, vre_ts, vre_params, bus_dict, vre_selection, vre_path):
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
        vre_dict[aid] = get_vre_dict(site, site_ts, vre_params, bus_dict, vre_selection)

    # write pickle
    out_file = vre_path
    utils.write_pickle(vre_dict, out_file)

def main(asset_path, ts_path, vre_selection, vre_path):
    '''
    This script takes creates the VRE dictionary for wind or solar assets for PyPSA_BC
    '''
    # Read in configuration file
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config2.yaml"
    cfg = utils.load_config(config_file)

    start_time = cfg['data']['cutout']['snapshots']['start'][0]
    end_time = cfg['data']['cutout']['snapshots']['end'][0]

    gen_generic = pd.read_csv(cfg['data']["coders"]["gen_generic"])
    vre_assets = pd.read_csv(asset_path)
    vre_ts = pd.read_csv(ts_path, index_col=0, parse_dates=True).loc[start_time:end_time]
    buses = pd.read_csv(cfg['output']['prepare_base_network']['folder'] + "/buses.csv")['name'].tolist()

    # (0A) Create folders if they have not been created already
    utils.create_folder(cfg['output']["pypsa_dict"]['folder'])

    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    bus_dict = utils.create_standard_gen_bus_map(buses)


    # (0C) Get vre cost information
    vre_params = get_vre_params(gen_generic, vre_selection)
                                
    # (1) Write pickle dictionaries for the vre assets.
    write_vre_dict(vre_assets, vre_ts, vre_params, bus_dict, vre_selection, vre_path)

    
if __name__ == '__main__':
    # main("results/interim/bc_wind_assets.csv", "results/interim/wind_ts.csv", "wind", "/home/pmcwhannel/repos/PyPSA_BC/results/pypsa-components/wind.pickle")
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])