import atlite
import numpy as np
import pandas as pd
import re
from pathlib import Path

'''
FUNCTIONS USED IN create_wind_assets.py
'''


def _infer_rated_kw_from_model(model: str):
    model = str(model)

    # Common form: MM114/3200 -> 3200 kW
    m = re.search(r"/(\d{3,5})\b", model)
    if m:
        return float(m.group(1))

    # Common form: V90/3MW -> 3000 kW
    m = re.search(r"/(\d+(?:\.\d+)?)\s*MW\b", model, flags=re.IGNORECASE)
    if m:
        return float(m.group(1)) * 1000.0

    # Common form: GE3.2-103 -> 3.2 MW
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*-", model)
    if m and float(m.group(1)) <= 30:
        return float(m.group(1)) * 1000.0

    return None


def _fallback_local_windturbine_config(model: str):
    """
    Fallback for turbine models not present in OEDB: choose the closest local
    atlite turbine by rated power (favor non-offshore names).
    """
    rated_kw = _infer_rated_kw_from_model(model)
    candidates = []
    for name, cfg_val in atlite.resource.windturbines.items():
        if isinstance(cfg_val, (str, Path)):
            cfg = atlite.resource.get_windturbineconfig(name)
        else:
            cfg = cfg_val
        p_kw = float(cfg["P"]) * 1000.0
        penalty = 5e5 if "offshore" in name.lower() else 0.0
        score = (abs(p_kw - rated_kw) if rated_kw is not None else 0.0) + penalty
        candidates.append((score, name, cfg))

    if not candidates:
        raise RuntimeError("No local atlite wind turbine configurations available for fallback.")

    _, name, cfg = min(candidates, key=lambda x: x[0])
    return cfg

#Get power capacity P for a wind turbine model, which is required for calculating the install capacity for a wind farm
#Used in generate_wind_assets()
#config_oedb = config_oedb in unique_models data frame
def get_power_cap(config_oedb):
    config = config_oedb.split('*')
    #ID used when there are multiple turbines with the same name, else just leave it blank
    # print(config[0])
    if len(config) > 1:
        if int(config[1]) >= 0:
            add = atlite.resource.get_oedb_windturbineconfig(config[0], id=int(config[1]))
            return add['P']

    
    try:
        add = atlite.resource.get_oedb_windturbineconfig(config[0])
    except RuntimeError as err:
        # Some model names map to multiple OEDB entries; fall back to first id.
        if "use `id` for an unambiguous search" in str(err):
            ids = [
                int(x)
                for x in re.findall(r"^\s*\d+\s+(\d+)\s+", str(err), flags=re.MULTILINE)
            ]
            if ids:
                add = atlite.resource.get_oedb_windturbineconfig(config[0], id=ids[0])
            else:
                raise
        elif "No turbine found" in str(err):
            add = _fallback_local_windturbine_config(config[0])
        else:
            raise
    return add['P']



'''
FUNCTIONS USED IN create_wind_ts.py
'''

#Function to return wind speed of closest pixel for a turbine
#Used in get_wind_coords()
#row = Some row in the wind_assets.csv data frame
#xaxis = Linear space ranging from westmost point on wind_atlas to eastmost point on wind_atlas
#yaxis = Linear space ranging from northmost point on wind_atlas to southmost point on wind_atlas
#data = The wind_atlas .tif data
def get_speed(row, xaxis, yaxis, data):
    #Get indices of the nearest pixels
    xIdx = np.searchsorted(xaxis, row['longitude'], side='left')
    yIdx = len(yaxis) - np.searchsorted(yaxis, row['latitude'], side='left', sorter=np.arange(len(yaxis)-1, -1, -1))

    return data[yIdx][xIdx] #Return the wind speed at the indices

#Generate a data frame that matches wind speeds from Global Wind Atlas to latitude/longitude values for scaling the cutout speeds
#Used in main()
#assets = The data frame for wind_assets.csv
#wind_atlas = The Global Wind Atlas wind speed data from the .tif file
#wind_geojson = The Global Wind Atlas geojson data which creates the shape for BC
def get_wind_coords(assets, wind_atlas, wind_geojson, province = 'MB'):
    # NOTE: Not consistently working across different provinces.. related to bounding box based on gwa_geojson. Current solution is the added "province" argument.

    # Store longitude and latitude values in a list for processing.
    # longitudes = [wind_geojson[i][0][j][0] for i in range(len(wind_geojson)) for j in range(len(wind_geojson[i][0]))] #[lon, lat], choose index 0
    # latitudes = [wind_geojson[i][0][j][1] for i in range(len(wind_geojson)) for j in range(len(wind_geojson[i][0]))] #[lon, lat], choose index 1
    if (province == 'MB') or (province == 'BC') :
        latitudes = []
        longitudes = []
        # 1st index: Changes latitude
        for row in wind_geojson:
            for coords in row[0]:
                longitudes.append(coords[0])
                latitudes.append(coords[1])

    else:
        longitudes = [wind_geojson[0][i][0] for i in range(len(wind_geojson[0]))]
        latitudes = [wind_geojson[0][j][1] for j in range(len(wind_geojson[0]))] # [lon, lat], choose index 1

    #Get latitude and longitude values to construct a bounding box for the wind speed data in latitude longitude format
    west = min(longitudes); north = max(latitudes) #Upper left corner
    east = max(longitudes); south = min(latitudes) #Lower right corner

    #Get x and y axis as linearly spaced longitudes and latitudes from the values calculated above
    xaxis = np.linspace(west, east, wind_atlas.shape[1])
    yaxis = np.linspace(north, south, wind_atlas.shape[0])

    #Match speeds of turbines to Global Wind Atlas
    wind_coords = assets.apply(lambda x: get_speed(x, xaxis, yaxis, wind_atlas), axis=1)

    return wind_coords


def get_XY(row, wnd):
    '''
    Function to get INDEX values of the square in the ERA5 data array is
    Used in generate_wind_ts()
    row = Some row in the wind_assets.csv data frame
    wnd = cutout.wnd100m.data
    '''
    x = 0
    y = 0
    for i in range(wnd.x.size):
        if row['x'] == wnd.x.values[i]:
            x = i
            break
    
    for j in range(wnd.y.size):
        if row['y'] == wnd.y.values[j]:
            y = j
            break

    return [x, y]


def scale_wind(row, wnd):
    '''
    Function to scale the wind speeds on the ERA5 data array
    Used in generate_wind_ts()
    row: Some row in the wind_assets.csv data frame
    wind: cutout.data.wnd100m
    NOTE: Modications made here 2023-10-25, since the flag parameter should not be used to dicate
          whether scaling occurs. Now the GWA scaling is used by default.
    '''

    wind_at_location = wnd.sel(x=row['x'], y=row['y']).values
    scaled = wind_at_location * row['GWA wind speed'] / np.mean(wind_at_location)
    return scaled

    
#Function for obtaining the wind turbine config from the wind_assets data frame
#Used in generate_wind_ts()
#config_oedb = config_oedb in unique_models data frame
#h = Hub height for a turbine farm
def get_config(config_oedb, h):
    # config = config_oedb.split('*')
    # #ID used when there are multiple turbines with the same name, else just leave it blank
    # if int(config[1]) >= 0:
    #     add = atlite.resource.get_oedb_windturbineconfig(config[0], id=int(config[1]))
    # else:
    #     add = atlite.resource.get_oedb_windturbineconfig(config[0])
    
    # add['hub_height'] = h #This hub height affects generation when cutout.wind

    config = config_oedb.split('*')
    #ID used when there are multiple turbines with the same name, else just leave it blank
    # print(config[0])
    if len(config) > 1:
        if int(config[1]) >= 0:
            add = atlite.resource.get_oedb_windturbineconfig(config[0], id=int(config[1]))
        else:
            try:
                add = atlite.resource.get_oedb_windturbineconfig(config[0])
            except RuntimeError as err:
                if "use `id` for an unambiguous search" not in str(err):
                    if "No turbine found" in str(err):
                        add = _fallback_local_windturbine_config(config[0])
                    else:
                        raise
                else:
                    ids = [
                        int(x)
                        for x in re.findall(r"^\s*\d+\s+(\d+)\s+", str(err), flags=re.MULTILINE)
                    ]
                    if not ids:
                        raise
                    add = atlite.resource.get_oedb_windturbineconfig(config[0], id=ids[0])
    else:
        try:
            add = atlite.resource.get_oedb_windturbineconfig(config[0])
        except RuntimeError as err:
            if "use `id` for an unambiguous search" in str(err):
                ids = [
                    int(x)
                    for x in re.findall(r"^\s*\d+\s+(\d+)\s+", str(err), flags=re.MULTILINE)
                ]
                if not ids:
                    raise
                add = atlite.resource.get_oedb_windturbineconfig(config[0], id=ids[0])
            elif "No turbine found" in str(err):
                add = _fallback_local_windturbine_config(config[0])
            else:
                raise

    add['hub_height'] = h #This hub height affects generation when cutout.wind
    return add