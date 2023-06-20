import yaml
import os
import atlite
import shapely
import geopandas as gpd
import pypsa
import numpy as np
import pandas as pd
import json
from shapely.ops import unary_union


def merge_assets(df,subset,sum_list):
    '''
    Function used to reduce hydroelectric datasets from turbines to an aggregate asset.
    This aggregation operator is currently applied only to the installed capacities for the units and
    the annual_avg_energy.
    df: Dataframe which is passed via a groupby operation.
    subset: Name of columns to use for deduplication
    sum_list: Name of parameters/columns to aggegtate using the sum operation.
    Example:
    Input dataframe has following entries below:
    component_id | asset_id | capacity | annual_avg_energy
    BC_MCA01_GEN | BC_MCA_GSS | 492 | 1936.79
    BC_MCA02_GEN | BC_MCA_GSS | 492 | 1936.79
    BC_MCA03_GEN | BC_MCA_GSS | 494 | 1942.7
    BC_MCA04_GEN | BC_MCA_GSS | 494 | 1942.7
    BC_MCA05_GEN | BC_MCA_GSS | 500 | 1968.45
    BC_MCA06_GEN | BC_MCA_GSS | 500 | 1968.45

    Output dataframe will have the following:
    asset_id | capacity | annual_avg_energy
    BC_MCA_GSS | 2972 | 11695.88

    '''
    # Other columns don't matter here for calcualting inflow and associated power production.
    df_out = df.drop_duplicates(subset=subset).set_index("connecting_node_code").copy()
    for param in sum_list:
        df_out[param] = df[param].sum()
    return df_out

def load_config(config_file):
    '''
    This function loads the configuration file for PyPSA_BC
    config_file: Path + filename of the configuration file. (i.e. ../config/config.yaml)
    '''
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)
    return cfg

def create_folder(folder):
    '''
    This functions creates a folder if not already created.
    If the folder is already created it takes no action
    folder: Path + folder name.
    '''
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Created folder @ {folder}")

def get_region_polygon(geometry):
    '''
    This function finds a bounding box of the region and creates a polygon for it.
    Returns a polygon of the regions max/min bounds in terms of lats and lons.
    '''
    if len(geometry) == 1:
        west_lon= geometry.bounds['minx'].iloc[0]
        south_lat  = geometry.bounds['miny'].iloc[0]
        east_lon = geometry.bounds['maxx'].iloc[0]
        north_lat = geometry.bounds['maxy'].iloc[0]
        bbox = (west_lon, south_lat, east_lon, north_lat)
        polygon = shapely.geometry.box(*bbox, ccw=True)
    else:
        print('There remains multiple geometries')
        exit(1)
    
    return polygon

def get_bounds(polygon_list):
    '''
    This function takes in a list of polygons and returns the maximum bounds for them.
    '''
    bounds = {}
    merged_poly = gpd.GeoSeries(unary_union(polygon_list))
    bounds["west_lon"] = merged_poly.geometry.bounds['minx'].iloc[0]
    bounds["south_lat"] = merged_poly.geometry.bounds['miny'].iloc[0]
    bounds["east_lon"] = merged_poly.geometry.bounds['maxx'].iloc[0]
    bounds["north_lat"] = merged_poly.geometry.bounds['maxy'].iloc[0]
    
    return bounds

def get_cutout_path(cfg):
    '''
    This function return the unique name based on the region and start/end year
    for a cutout. 
    return: file path + name for the cutout described by selections in the
    cutout configuration.
    '''
    # Create file_path name with custom year_date
    start_year = cfg['cutout']["snapshots"]["start"][0][:4]
    end_year = cfg['cutout']["snapshots"]["end"][0][:4]
    prefix = cfg['cutout']['path'] + cfg['cutout']['region']["name"]

    if start_year == end_year:
        suffix = start_year
        file = "_".join([prefix, suffix + ".nc"])
    else: # multi_year_file
        suffix = "_".join([start_year, end_year])
        file = "_".join([prefix, suffix + ".nc"])

    return file

def create_era5_cutout(bounds, cfg):
    '''
    This function creates a cutout based on data for era5.
    '''
    # Extract parameters from configuration file
    dx,dy = cfg["cutout"]["dx"], cfg["cutout"]['dy']
    time_horizon = slice(cfg["cutout"]["snapshots"]['start'][0],
                        cfg["cutout"]["snapshots"]['end'][0])

    # get path + filename for the cutout
    file = get_cutout_path(cfg)

    # Create the cutout based on bounds found from above
    cutout = atlite.Cutout(path=file,
                    module=cfg["cutout"]["source"],
                    x=slice(bounds['west_lon'] - dx, bounds['east_lon'] + dx),
                    y=slice(bounds['south_lat'] - dy, bounds['north_lat'] + dy ),
                    dx=dx,
                    dy=dy,
                    time=time_horizon)

    cutout.prepare()

def convert_cid_2_aid(cid,old_aid):
    '''
    This creates an asset id (aid) based on the component id (cid).
    Common Example: 
            cid -> BC_ZBL03_GEN
            old_aid -> BC_ZBL_GSS
            new_aid -> BC_ZBL_GSS
           Example:
           cid -> BC_BR0101_GEN
           old_aid -> BC_BR1_DFS
           new_aid -> BC_BR1_DFS
           Example:
           cid -> BC_BR0102_GEN
           old_aid -> BC_BR2_GSS
           new_aid -> BC_BR2_GSS

    '''
    aid_start = old_aid.split('_')[0]
    cid_2_aid = cid.split('_')[1][:3]
    aid_end = old_aid.split('_')[-1]
    
    if aid_start != cid.split('_')[0]: # error check
        print('Error detected in convert_cid_2_aid')
        exit(3)
    new_aid= "_".join([aid_start, cid_2_aid, aid_end])
    return new_aid 

def write_pickle(data_dict, filepath):
    '''
    Write a pickle file based on a dictionary.
    '''
    with open(filepath,"wb") as f:
        pd.to_pickle(data_dict, f)
    f.close()
    print(f'Wrote pickle file {filepath}')

def read_pickle(filepath):
    '''
    Read a json file based on a dictionary.
    '''
    with open(filepath, 'rb') as f:
        data_dict = pd.read_pickle(f) 
    return data_dict

def create_standard_gen_bus_map(buses):
    '''
    This function accepts a list a buses and returns a mapping from the buses to the lowest voltage
    bus for each unique node code. The underlying assumption used in selecting a bus to connect a generator to is that
    generators are connected the lowest voltage bus at their given node location.
    Example:
    Buses = ["230_ABN_GSS", "138_ABN_GSS","500_MCA_GSS", "63_MCA_GSS"] -> bus_dict = {ABN_GSS:138, MCA_GSS:63}
    '''
    bus_dict = {}
    # Assume generators are connected to lowest voltage bus at their given node_code
    for bus in buses:
        node = "_".join(bus.split('_')[1:])
        voltage = int(bus.split('_')[0])
        if node not in bus_dict.keys():
            bus_dict[node] = voltage
        else:
            bus_dict[node] = min(voltage,bus_dict[node])

    return bus_dict

def get_gen_bus(node_code, bus_dict):
    '''
    This function returns the correct standardized electric bus for a generator.
    Example:
    node_code = "BC_ABS_GSS" 
    bus_dict = {"ABS_GSS":230, "MCA_GSS":63}
    return -> "230_ABS_GSS"
    '''
    node_code_suffix = "_".join(node_code.split('_')[1:]) # i.e.  BC_ABN_GSS -> ABN_GSS
    return "_".join([str(bus_dict[node_code_suffix]), node_code_suffix])

def get_multi_link_override():
    '''
    Gets the multi-link override. Needed for cascaded hydroelectric.
    '''
    # From PyPSA CHP Example: This ensures we can add 2 outputs for a single link i.e bus0 -> bus_1 AND bus_2
    override_component_attrs = pypsa.descriptors.Dict(
        {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    )
    override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "2nd bus",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        1.0,
        "2nd bus efficiency",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "2nd bus output",
        "Output",
    ]
    return override_component_attrs