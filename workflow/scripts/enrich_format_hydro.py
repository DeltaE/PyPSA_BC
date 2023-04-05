import pypsa
from bc_power import utils
import pandas as pd

def get_bus_dict(bus_name):
    '''
    This function is used to create a dictionary for a water bus of a cascade or hydro reservoir.
    The water bus is used as medium to move water through the reservoir for
    discharge and power production, spillage into another water bus, or to the reservoirs storage.
    OR
    This function is used to create a reservoir bus for cascade.
    The reservoir bus is connected to connect to a store in pypsa which is used
    to represent the storage of a water reservoir for large controllable dams (mainly hydroelectric dams).
    '''
    bus_dict = {"class_name":"Bus",
                "name":bus_name,
                "carrier":"water"
                }
    
    return bus_dict

def get_res_store_dict(link_name, bus0_name, bus1_name, p_nom):
    '''
    This function is used to get a reservoir storage link dictionary.
    '''

    store_dict =  {"class_name":"Link",
                    "name": link_name,
                    "bus0": bus0_name,
                    "bus1": bus1_name,
                    "efficiency":1., # mass balance
                    "p_nom":1e100, # UPDATE: Should be derived to ensure larger than max(inflow + spill + discharge of any upstream reservoirs)
                    }
    
    return store_dict


def get_res_release_dict(link_name, bus0_name, bus1_name, p_nom):
    '''
    This function is used to get a reservoir release link dictionary.
    p_nom: Should be equal to storage capacity + discharge + spill (adjusted for timesteps...)
    '''

    store_dict =  {"class_name":"Link",
                    "name": link_name,
                    "bus0": bus0_name,
                    "bus1": bus1_name,
                    "efficiency":1., # mass balance
                    "p_nom":p_nom, # max possible release from reservoir to water bus.
                    }
    
    return store_dict

def main():
    '''
    This script takes the hydro generation files and formats them into dataframes which will be used
    to add network components in PyPSA along with joined by their inflow or RoR power availability
    counterparts.
    The inputs to this script are the following:
    1) Hydro Generation
    2) Hydro Reservoirs

    The output of this script are the following 2 dataframes:
    1) A dataframe of hydro assets to be modelled as RoR assets with their corresponing attributes.
    2) A dataframe of hydro assets to be modelled as reservoir assets with their corresponding attributes.
    '''
    # Determine
    # read config for paths
    # Read in configuration file
    config_file = r"/mnt/c/Users/pmcw9/Delta-E/PICS/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    hydro_gen = pd.read_csv(cfg["hydro_prep"]["hydro_generation"])
    hydro_res = pd.read_csv(cfg["hydro_prep"]["hydro_reservoir"])


if __name__ == '__main__':
    main()