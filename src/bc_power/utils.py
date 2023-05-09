import yaml
import os
import atlite
import shapely
import geopandas as gpd
import pandas as pd
import json
from shapely.ops import unary_union

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

def create_era5_cutout(bounds, cfg):
    '''
    This function creates a cutout based on data for era5.
    '''
    # Extract parameters from configuration file
    dx,dy = cfg["cutout"]["dx"], cfg["cutout"]['dy']
    time_horizon = slice(cfg["cutout"]["snapshots"]['start'][0],
                        cfg["cutout"]["snapshots"]['end'][0])

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
