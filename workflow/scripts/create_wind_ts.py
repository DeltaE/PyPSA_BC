import sys
from pathlib import Path
import rasterio as rio
import pandas as pd
import numpy as np
import geojson as gj
import atlite
import geopandas as gpd
import xarray as xr

'''
====================================================================================================

    HELPER FUNCTIONS BELOW HERE, TO BE MOVED TO A MODULE EVENTUALLY

====================================================================================================
'''

'''
    get_wind_coords() SECTION
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



'''
    generate_wind_ts() SECTION
'''

#Function to get INDEX values of the square in the ERA5 data array is
#Used in generate_wind_ts()
#row = Some row in the wind_assets.csv data frame
#wind = cutout.wnd100m.data
def get_XY(row, wind):
    x = 0
    y = 0
    for i in range(wind.x.size):
        if row['x'] == wind.x.values[i]:
            x = i
            break
    
    for j in range(wind.y.size):
        if row['y'] == wind.y.values[j]:
            y = j
            break

    return [x, y]

#Function to scale the wind speeds on the ERA5 data array
#Used in generate_wind_ts()
#row = Some row in the wind_assets.csv data frame
#wind = cutout.wnd100m.data
def scale_wind(row, wind):
    if row['Flag'] == 1:
        #Scale the wind speeds at this location on the ERA5 data array
        wind_at_location = wind.sel(x=row['x'], y=row['y']).values
        scaled = wind_at_location * row['GWA wind speed'] / np.mean(wind_at_location)
        return scaled
    else:
        #Do nothing
        return None

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


#Function to calculate the wind generation profiles from the cutout, for a given turbine farm
#Used in generate_wind_ts()
#cutout = Cutout object
#cluster = Turbine cluster, these are rows from wind_assets
def calculate_wind(cutout, cluster):
    cap_factors = cutout.wind(turbine=get_config(cluster['config_oedb'], cluster['Hub height (m)']), capacity_factor=True)

    # This code calculates the cells/grids for where generation exists
    # Cells of the cutout
    cells = cutout.grid 

    #site: location of the turbines with latitudes and longitudes specified
    sites = gpd.GeoDataFrame([[cluster['location'], cluster['longitude'], cluster['latitude'], cluster['Install capacity']]],
                            columns=['name', 'lon', 'lat', 'capacity']
                            ).set_index('name')

    # Finds cutout cells nearest to (lat,lon) of each site 
    # (x,y) of cells/grids are the center points by default and correspong to (lat,lon) of the cells.
    nearest = cutout.data.sel(
        {'x': sites.lon.values, 'y': sites.lat.values}, 'nearest').coords

    # Add new columns for the grid (x,y) where each generator falls.
    sites['x'] = nearest.get('x').values
    sites['y'] = nearest.get('y').values

    # Performs an inner join on the (x,y) values
    cells_generation = sites.merge(
        cells, how='inner').rename(pd.Series(sites.index))

    # Layout indicates the installed capacity and grid location of installed capacity
    layout = xr.DataArray(cells_generation.set_index(['y', 'x']).capacity.unstack())\
                        .reindex_like(cap_factors).rename('Installed Capacity [MW]')

    # Finally generate the power produced at these sites
    # Returns the power generation of the generators for each location
    # The data is returned as an xarray with shape of (timesteps, generators)
    power_generation = cutout.wind(turbine=get_config(cluster['config_oedb'], cluster['Hub height (m)']),
                                    layout = layout,
                                    shapes = cells_generation.geometry,
                                    per_unit = False) # currently using per-unit
    
    return power_generation.to_pandas()


'''
    calibrate_generation() SECTION
'''

#Function to calculate some value 'a' such that sum(a*w) = P_rated
#Used in calibrate_generation()
#w_vec = Vector of wind generation values
#aag = Target AAG to calibrate to
#p = P_rated
#a_init = Initial guess for the value of 'a'
def getA(w_vec, aag, p, a_init):
    #Initialize values for the loop
    iter = 0
    a = a_init

    while iter < 100: #Stop after 100 iterations if it does not work out for some reason (this should never happen in practice)
        #dfda:
        #-w, when 0 <= a*w <= P_rated
        #0, everywhere else
        dfda = w_vec.apply(lambda w: float(a*w >= 0) * float(a*w <= p) * w * -1).sum()

        #The value for f = AAG - sum of t_{i} where t_{i} = min(a*w_i, P_rated)
        f = aag - w_vec.apply(lambda w: min(a*w, p)).sum()

        #Get new value for 'a', a_{i+1} = a_{i} - (f / dfda)
        a_new = a - (f / dfda)

        #Calculate tolerance percentage
        tol = abs(100 * ((a_new - a) / a))

        #What to do based on the value of the tolerance
        if tol < 0.01:
            break #Exit loop and return aNew if tolerance is small enough
        else:
            a = a_new #Loop again using the newly calculated 'a' value

        iter += 1
    
    return a_new

#Function to calibrate the wind generation data by using the CODERS AAG data
#Used in main()
#wind_ts = Wind time series data frame from generate_wind_ts()
#wind_assets = Wind assets data frame from wind_assets.py, this has the CODERS AAG data ready
def calibrate_generation(wind_ts, wind_assets):
    #Get date breakpoints for each year in the wind time series
    years = wind_ts.groupby(wind_ts.index.year).count()[wind_ts.columns[0]]
    year_index = [0]

    #Build year_index, which will be used to calibrate wind generations per year
    for year in years.index:
        year_index.append(year_index[-1] + years.loc[year])
    
    #Empty DataFrame whose columns will be filled
    calibrated_wind = pd.DataFrame(data={'time': wind_ts.index})

    #Generate calibrated wind generation time series
    for component in wind_assets['component_id'].unique():
        #Get necessary parameters for calculating the value 'a' for calibration
        wind = wind_ts[component].reset_index()
        aag = wind_assets.loc[wind_assets['component_id'] == component]['CODERS AAG (GWh/y)'].values[0] * 1000 #GWh --> MWh
        p = wind_assets.loc[wind_assets['component_id'] == component]['Install capacity'].sum()
        
        #Empty series to hold the wind generation for a component_id
        column_to_add = pd.Series(name=component, dtype='float64')

        #Calculate calibrated wind generations by year
        for i in range(len(year_index) - 1):
            subset = wind.loc[(wind.index >= year_index[i]) & (wind.index < year_index[i+1])].copy()

            #Obtain a value of a for the year
            a = getA(subset[component], aag, p, aag / subset[component].sum())

            #Use the a value to normalize the wind generation in that year
            subset[component] = subset[component].apply(lambda w: min(a*w, p))

            column_to_add = pd.concat([column_to_add, subset[component]], ignore_index=True)
        
        #Add the column to the final data drame
        calibrated_wind = pd.concat([calibrated_wind, column_to_add], axis=1)
    

    #Final calibrated wind generation data frame to return, set index to time to match the uncalibrated data frame's format
    return calibrated_wind.set_index('time')

'''
====================================================================================================

    THE ACTUAL FUNCTIONS THAT PRODUCE wind_ts.csv ARE BELOW HERE

====================================================================================================
'''

#This part builds the wind_ts data frame to be written into a CSV file (NO CALIBRATION WITH CODERS AAG DATA)
#Used in main()
#wind_assets = Wind assets data frame from wind_assets.csv, after a bit of modification from get_wind_coords()
#cutout_path = Path to the cutout file
def generate_wind_ts(wind_assets, cutout_path):
    #Load in the cutout for scaling with Global Wind Atlas wind speeds stored in wind_assets
    cutout = atlite.Cutout(path=cutout_path)

    #The wind data is the main target here
    wind = cutout.data.wnd100m


    #Select the nearest squares on the ERA5 grid for each wind turbine farm
    nearest = cutout.data.sel({'x': wind_assets.longitude.values, 'y': wind_assets.latitude.values}, 'nearest').coords
    x_near = nearest.get('x').values
    y_near = nearest.get('y').values

    #Put these matched squares into the dataframe to line them up nice
    wind_assets['x'] = x_near
    wind_assets['y'] = y_near

    #Get scaled wind values and x-y coords
    scaled_wind = wind_assets.apply(lambda x: scale_wind(x, wind), axis=1)
    xy = wind_assets.apply(lambda x: get_XY(x, wind), axis=1)

    #Combine scaled_wind and xy an reformat a bit
    scaled_wind_xy = pd.DataFrame(data={'Wind': scaled_wind, 'XY': xy}).dropna().reset_index().drop(columns='index')
    
    #This will be used to hold the data array values to be assigned back to the cutout's dataset
    to_data_array = wind

    #Writing the data in scaledWindXY into toDataArray
    for i in range(scaled_wind_xy.index.size):
        #to_data_array[time, y, x], scaledWind.Wind contains the scaled hourly wind at location y, x
        to_data_array[:, scaled_wind_xy.loc[i].XY[1], scaled_wind_xy.loc[i].XY[0]] = scaled_wind_xy.loc[i].Wind

    #Overwriting the wind data in the cutout with the scaled data
    wind.data = to_data_array


    #Generating the wind generation time series here
    #Loop through wind_assets and construct the wind_generation dataframe
    wind_generation = pd.DataFrame()

    for i in range(wind_assets.index.size):
        wind_generation[str(wind_assets.index[i])] = calculate_wind(cutout, wind_assets.loc[wind_assets.index[i]])
    
    #Rename columns to their respective asset_id
    wind_generation.columns = wind_assets['asset_id'].values


    #Summing values for any component_id that has more than one asset_id associated with it
    for i in range(wind_assets.index.size):
        #A 'Flag' value of 0 means that this asset_id is part of a component_id
        if wind_assets.loc[i]['Flag'] == 0:
            #The main column that will store the value for the sum
            column_to_add_to = wind_assets.loc[i]['component_id'] + '01'

            #Compute the sum to the main column, then drop the column of the asset_id
            wind_generation[column_to_add_to] = wind_generation.apply(lambda x: x[column_to_add_to] + x[wind_assets.loc[i]['asset_id']], axis=1)
            wind_generation = wind_generation.drop(columns=[wind_assets.loc[i]['asset_id']])

    #Finish by renaming the columns to their respective component_id names
    wind_generation.columns = wind_assets['component_id'].unique()


    return wind_generation

#Does some input verification before generating the time series
def main():
    try:
        #Try reading the arguments passed in the terminal
        assets_path = Path(sys.argv[1])
        cutout_path = Path(sys.argv[2])
        wind_atlas_path = Path(sys.argv[3])
        wind_geojson_path = Path(sys.argv[4])
        calibration_flag = int(sys.argv[5]) #0 for no calibration, 1 for calibration
        output_path = Path(sys.argv[6])
    except:
        #Less than 6 arguments, return error code 1
        print('There are inputs missing')
        return 1
    else:
        #Correct number of arguments
        try:
            #Load the wind_assets
            assets = pd.read_csv(assets_path)

            #Load in the Global Wind Atlas wind speeds for BC with rasterio
            with rio.open(wind_atlas_path) as f:
                wind_atlas = f.read(1)
                f.close()

            #Load in the Global Wind Atlas geojson
            with open(wind_geojson_path) as f:
                wind_geojson = gj.load(f)['geometry']['coordinates']
                f.close()

        except:
            #An input may be spelled incorrectly, return error code 2
            print('One or more inputs are in the wrong format')
            return 2
        else:
            #Start by appending the Global Wind Atlas wind speeds to assets
            assets['GWA wind speed'] = get_wind_coords(assets, wind_atlas, wind_geojson)

            #All is good, start generating wind_ts data frame
            if calibration_flag == 0:
                #No calibration for the wind speeds with CODERS AAG
                wind_ts = generate_wind_ts(assets, cutout_path)

            elif calibration_flag == 1:
                #Use CODERS AAG to calibrate the wind speeds
                wind_ts = calibrate_generation(generate_wind_ts(assets, cutout_path), assets)

            else:
                #Calibration flag, return error code 3
                print('Invalid value for calibration flag, choose 0 for no calibration or 1 for calibration')
                return 3
            
            #Write wind_ts.csv
            wind_ts.to_csv(output_path)

            #Return code 0 is for when everything runs without a problem
            return 0

if __name__ == '__main__':
    main()
    