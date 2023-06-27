import sys
from pathlib import Path
import pandas as pd
import atlite
from bc_power import solar_wind

#This part builds the solar_ts data frame to be written into a CSV file (NO CALIBRATION WITH CODERS AAG DATA)
#Used in main()
#solar_assets = Solar assets data frame from solar_assets.csv
#cutout_path = Path to the cutout file
def generate_solar_ts(solar_assets, cutout_path):
    #Load in the cutout
    cutout = atlite.Cutout(path=cutout_path)


    #Generating the solar generation time series here
    #Loop through solar_assets and construct the solar_generation dataframe
    solar_generation = pd.DataFrame()

    for i in range(solar_assets.index.size):
        solar_generation[str(solar_assets.index[i])] = solar_wind.calculate_MW(cutout, solar_assets.loc[solar_assets.index[i]], 'solar')
    
    #Rename columns to their respective asset_id
    solar_generation.columns = solar_assets['asset_id'].values


    #Summing values for any component_id that has more than one asset_id associated with it
    for i in range(solar_assets.index.size):
        #A 'Flag' value of 0 means that this asset_id is part of a component_id
        if solar_assets.loc[i]['Flag'] == 0:
            #The main column that will store the value for the sum
            column_to_add_to = solar_assets.loc[i]['component_id'] + '01'

            #Compute the sum to the main column, then drop the column of the asset_id
            solar_generation[column_to_add_to] = solar_generation.apply(lambda x: x[column_to_add_to] + x[solar_assets.loc[i]['asset_id']], axis=1)
            solar_generation = solar_generation.drop(columns=[solar_assets.loc[i]['asset_id']])

    #Finish by renaming the columns to their respective component_id names
    solar_generation.columns = solar_assets['component_id'].unique()


    return solar_generation

#Does some input verification before generating the time series
def main():
    try:
        #Try reading the arguments passed in the terminal
        assets_path = Path(sys.argv[1])
        cutout_path = Path(sys.argv[2])
        calibration_flag = int(sys.argv[3]) #0 for no calibration, 1 for calibration
        output_path = Path(sys.argv[4])
    except:
        #Less than 4 arguments, return error code 1
        print('There are inputs missing')
        return 1
    else:
        #Correct number of arguments
        try:
            #Load the wind_assets
            assets = pd.read_csv(assets_path)

        except:
            #An input may be spelled incorrectly, return error code 2
            print('One or more inputs are in the wrong format')
            return 2
        else:
            #All is good, start generating wind_ts data frame
            if calibration_flag == 0:
                #No calibration for the solar generation with CODERS AAG
                solar_ts = generate_solar_ts(assets, cutout_path)

            elif calibration_flag == 1:
                #Use CODERS AAG to calibrate the wind speeds
                solar_ts = solar_wind.calibrate_generation(generate_solar_ts(assets, cutout_path), assets)

            else:
                #Calibration flag, return error code 3
                print('Invalid value for calibration flag, choose 0 for no calibration or 1 for calibration')
                return 3
            
            #Write solar_ts.csv
            solar_ts.to_csv(output_path)

            #Return code 0 is for when everything runs without a problem
            return 0

if __name__ == '__main__':
    main()
    