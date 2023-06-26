import pypsa
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import pyomo
from bc_power import utils


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

def get_busmap_dict(network, gadm_bc):
    # (1) Determine mapping between buses and GADM_region
    # Determine which buses line within which gadm regions
    # output: {bus:GADM_region, ...} 
    busmap_dict = {}
    for bus,row in network.buses.iterrows(): # Loop over buses
        point = Point((row['x'],row['y']))
        for idx, row in gadm_bc.iterrows():
            if row.geometry.contains(point): # Check to see which GADM region is in contained within
                busmap_dict[bus] = row["NAME_2"]
                break

    return busmap_dict

def create_adm_buses(network, gadm_bc, busmap_dict):
    '''
    Created buses for the administrative regions.
    '''
    # (2) Create a new bus for GADM regions which have been mapped
    for gadm_region in busmap_dict.values(): 
        if gadm_region not in network.buses.index.to_list():
            # Could also use a dictionary to provide these values
            network.add(class_name="Bus",
                        name=gadm_region,
                        x=gadm_bc[gadm_bc["NAME_2"] == gadm_region].geometry.centroid.x.iloc[0],
                        y=gadm_bc[gadm_bc["NAME_2"] == gadm_region].geometry.centroid.y.iloc[0],
                        v_nom = 300
                        )
        else:
            pass

def replace_bus_refs(component, col, busmap_dict):
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
        network.remove(class_name='Bus',name=bus_name)


    # (iii) remove transformers
    trans_to_remove = []
    for trans_name,row in network.transformers.iterrows():
        trans_to_remove.append(trans_name)

    for trans_name in trans_to_remove:
        network.remove(class_name='Transformer',name=trans_name)

    # NOTE: Custom revmoval of Site C
    # Will need to be updated
    network.remove(class_name='Link', name='BC_STC_GSS Discharge Link')



def main():
    '''
    This script is used to build the model.(Currently, designed to build the existing electricity system in BC. (with site-c))
    The scripts takes in the followig data:
    
    1) Network structure
    2) Hydro assets
    3) Hydro timeseries data
    4) Wind assets
    5) Wind timeseries data
    6) Solar assets
    7) Solar timeseries data

    The script will save a NetCDF file of the instantiated network. 

    '''
    # (0) Load config file
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    # (1) Load files
    network = pypsa.Network(override_component_attrs=utils.get_multi_link_override())
    network_path = cfg['network']['folder'] 
    hydro_ror_path = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"]["hydro_ror"]
    hydro_res_path = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"]["hydro_res"]
    wind_path = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"]["wind"]
    pv_path = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"]["solar"]
    tpp_path = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"]["tpp"]

    bus_path = cfg['network']['folder'] + "/buses.csv" # buses.csv should never be changed 


    network.import_from_csv_folder(network_path)
    ror_dict = utils.read_pickle(hydro_ror_path)
    res_dict = utils.read_pickle(hydro_res_path)
    wind_dict = utils.read_pickle(wind_path)
    pv_dict = utils.read_pickle(pv_path)
    tpp_dict = utils.read_pickle(tpp_path)

    # (2) Set timeslicing
    network.set_snapshots(ror_dict['BC_ABN_GSS']['p_max_pu'].index) # UDPATE REQUIRED

    # (3) Start adding assets
    add_hydro_ror_assets(network, ror_dict)
    add_hydro_res_assets(network, res_dict)
    add_pv_assets(network, pv_dict)
    add_wind_assets(network, wind_dict)
    add_tpp_assets(network, tpp_dict)

    # (4) add carriers outside of default ELC
    network.add("Carrier","NG", co2_emissions=1.0)

    # (5) Clustering
    # NOTE: Centroid calculation will need an update and validation.
    
    bus_dict = {name:0 for name in pd.read_csv(bus_path)['name'].tolist()}
    geojson_file = "/mnt/c/Users/pmcw9/Delta-E/PICS/Data/regions/gadm41_CAN_2.json"
    gdf = gpd.read_file(geojson_file)

    # Get GeoDataFrame of the GADM regions.
    gadm_bc = gdf[gdf["NAME_1"] =="BritishColumbia"]

    # (6) get busmap dictionary
    busmap_dict = get_busmap_dict(network, gadm_bc)

    # (7) Add new gadm regions as buses
    create_adm_buses(network, gadm_bc, busmap_dict)

    # (8) Determine all componenets to relink
    replace_bus_refs(network.generators,'bus', busmap_dict)
    replace_bus_refs(network.links,'bus0', busmap_dict)
    replace_bus_refs(network.links,'bus1', busmap_dict)
    replace_bus_refs(network.links,'bus2', busmap_dict)
    replace_bus_refs(network.lines,'bus0', busmap_dict)
    replace_bus_refs(network.lines,'bus1', busmap_dict)
    replace_bus_refs(network.loads,'bus', busmap_dict)

    # (9) Remove the old components which are no longer needed (lines, buses)
    remove_old_components(network, busmap_dict)

    # NOTE: Think about how to aggregate the new lines...
    network.lines['s_nom'] = 15000 # Network (inf capacity basically)

    # (10) Add load

    start_time = cfg['scope']['temporal']['start'] 
    end_time = cfg['scope']['temporal']['end']

    res_load = pd.read_csv(r"/home/pmcwhannel/repos/PyPSA_BC/results/hourly_res.csv",
                            index_col=0, parse_dates=True).loc[start_time:end_time]
    csmi_load = pd.read_csv(r"/home/pmcwhannel/repos/PyPSA_BC/results/hourly_csmi.csv",
                            index_col=0, parse_dates=True).loc[start_time:end_time]
    
    # NOTE: Update coming for Stikine and CentralCoast.
    for col in res_load.columns: 
        if col not in ["Stikine", "CentralCoast", "NorthernRockies"]: 
            load_ts= res_load[col] + csmi_load[col]
            network.add("Load", "{} ELC Load".format(col), bus=col, p_set=load_ts)
        else:
            pass
    
    # Solve network
    network.optimize(solver_name='cbc')

    # Save network
    network.export_to_netcdf(cfg["pypsa"]["results"])
if __name__ == '__main__':
    main()