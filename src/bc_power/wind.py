import atlite
import numpy as np
import pandas as pd

'''
FUNCTIONS USED IN create_wind_assets.py
'''

#Get power capacity P for a wind turbine model, which is required for calculating the install capacity for a wind farm
#Used in generate_wind_assets()
#config_oedb = config_oedb in unique_models data frame
def get_power_cap(config_oedb):
    config = config_oedb.split('*')
    #ID used when there are multiple turbines with the same name, else just leave it blank
    if int(config[1]) >= 0:
        add = atlite.resource.get_oedb_windturbineconfig(config[0], id=int(config[1]))
    else:
        add = atlite.resource.get_oedb_windturbineconfig(config[0])

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
def get_wind_coords(assets, wind_atlas, wind_geojson):
    #Store longitude and latitude values in a list for processing.
    longitudes = [wind_geojson[i][0][j][0] for i in range(len(wind_geojson)) for j in range(len(wind_geojson[i][0]))] #[lon, lat], choose index 0
    latitudes = [wind_geojson[i][0][j][1] for i in range(len(wind_geojson)) for j in range(len(wind_geojson[i][0]))] #[lon, lat], choose index 1

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
    config = config_oedb.split('*')
    #ID used when there are multiple turbines with the same name, else just leave it blank
    if int(config[1]) >= 0:
        add = atlite.resource.get_oedb_windturbineconfig(config[0], id=int(config[1]))
    else:
        add = atlite.resource.get_oedb_windturbineconfig(config[0])
    
    add['hub_height'] = h #This hub height affects generation when cutout.wind

    return add