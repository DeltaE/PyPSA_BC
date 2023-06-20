import sys
from pathlib import Path
import rasterio as rio
import pandas as pd
import geojson as gj
import atlite
from bc_power import wind, solar_wind, utils

#This part builds the wind_ts data frame to be written into a CSV file (NO CALIBRATION WITH CODERS AAG DATA)
#Used in main()
#wind_assets = Wind assets data frame from wind_assets.csv, after a bit of modification from get_wind_coords()
#cutout_path = Path to the cutout file
def generate_wind_ts(wind_assets, cutout_path):
    #Load in the cutout for scaling with Global Wind Atlas wind speeds stored in wind_assets
    cutout = atlite.Cutout(path=cutout_path)

    #The wind data is the main target here
    wnd = cutout.data.wnd100m

    #Select the nearest squares on the ERA5 grid for each wind turbine farm
    nearest = cutout.data.sel({'x': wind_assets.longitude.values, 'y': wind_assets.latitude.values}, 'nearest').coords
    x_near = nearest.get('x').values
    y_near = nearest.get('y').values

    #Put these matched squares into the dataframe to line them up nice
    wind_assets['x'] = x_near
    wind_assets['y'] = y_near

    #Get scaled wind values and x-y coords
    scaled_wind = wind_assets.apply(lambda x: wind.scale_wind(x, wnd), axis=1)
    xy = wind_assets.apply(lambda x: wind.get_XY(x, wnd), axis=1)

    #Combine scaled_wind and xy an reformat a bit
    scaled_wind_xy = pd.DataFrame(data={'Wind': scaled_wind, 'XY': xy}).dropna().reset_index().drop(columns='index')
    
    #This will be used to hold the data array values to be assigned back to the cutout's dataset
    to_data_array = wnd

    #Writing the data in scaledWindXY into toDataArray
    for i in range(scaled_wind_xy.index.size):
        #to_data_array[time, y, x], scaledWind.Wind contains the scaled hourly wind at location y, x
        to_data_array[:, scaled_wind_xy.loc[i].XY[1], scaled_wind_xy.loc[i].XY[0]] = scaled_wind_xy.loc[i].Wind

    #Overwriting the wind data in the cutout with the scaled data
    wnd.data = to_data_array


    #Generating the wind generation time series here
    #Loop through wind_assets and construct the wind_generation dataframe
    
    wind_gen_dict = {} # new

    # for i in range(wind_assets.index.size):
    #     # Right here create a dictionary with key as asset_id and pd
    #     wind_generation[str(wind_assets.index[i])] = solar_wind.calculate_MW(cutout, wind_assets.loc[wind_assets.index[i]], 'wind')

    for _,row in wind_assets.iterrows():
        if row['asset_id'] not in wind_gen_dict.keys(): 
            wind_gen_dict[row['asset_id']] = solar_wind.calculate_MW(cutout, row, 'wind').squeeze()
        else:
            wind_gen_dict[row['asset_id']] += solar_wind.calculate_MW(cutout, row, 'wind').squeeze()

    wind_generation = pd.DataFrame(wind_gen_dict)
    # Rename columns to their respective asset_id
    # wind_generation.columns = wind_assets['asset_id'].values # should be a non-unique issue occuring here 


    # Summing values for any component_id that has more than one asset_id associated with it
    # for i in range(wind_assets.index.size):
    #     # A 'Flag' value of 0 means that this asset_id is part of a component_id
    #     if wind_assets.loc[i]['Flag'] == 0:
    #         #The main column that will store the value for the sum
    #         column_to_add_to = wind_assets.loc[i]['component_id'] + '01'

    #         #Compute the sum to the main column, then drop the column of the asset_id
    #         wind_generation[column_to_add_to] = wind_generation.apply(lambda x: x[column_to_add_to] + x[wind_assets.loc[i]['asset_id']], axis=1)
    #         wind_generation = wind_generation.drop(columns=[wind_assets.loc[i]['asset_id']])

    #Finish by renaming the columns to their respective connecting_node_code names
    # wind_generation.columns = wind_assets['connecting_node_code'].unique()


    return wind_generation

#Does some input verification before generating the time series
def main():
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    #Try reading the arguments passed in the terminal
    assets_path = cfg['wind']['asset_path'] # Path(sys.argv[1])
    cutout_path = utils.get_cutout_path(cfg) # Path(sys.argv[2])
    wind_atlas_path = cfg['wind']['atlas_speed'] # Path(sys.argv[3])
    wind_geojson_path = cfg['wind']['atlas_geojson'] # Path(sys.argv[4])
    calibration_flag = cfg['wind']['calibration'] # 0 for no calibration, 1 for calibration
    output_path = cfg['wind']['ts_path'] # Path(sys.argv[6])


    #Correct number of arguments
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


    #Start by appending the Global Wind Atlas wind speeds to assets
    assets['GWA wind speed'] = wind.get_wind_coords(assets, wind_atlas, wind_geojson)

    #All is good, start generating wind_ts data frame
    if calibration_flag == 0:
        #No calibration for the wind generation with CODERS Annual Average Generation
        wind_ts = generate_wind_ts(assets, cutout_path)

    elif calibration_flag == 1:
        #Use CODERS AAG to calibrate the wind generation
        wind_ts = solar_wind.calibrate_generation(generate_wind_ts(assets, cutout_path), assets)

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
    