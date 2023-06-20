import sys
from pathlib import Path
import pandas as pd
import atlite
from bc_power import solar_wind, utils

#This part builds the solar_ts data frame to be written into a CSV file (NO CALIBRATION WITH CODERS AAG DATA)
#Used in main()
#solar_assets = Solar assets data frame from solar_assets.csv
#cutout_path = Path to the cutout file
def generate_solar_ts(solar_assets, cutout_path):
    #Load in the cutout
    cutout = atlite.Cutout(path=cutout_path)

    #Generating the solar generation time series here
    #Loop through solar_assets and construct the solar_generation dataframe
    wind_gen_dict = {} # new

    # for i in range(wind_assets.index.size):
    #     # Right here create a dictionary with key as asset_id and pd
    #     wind_generation[str(wind_assets.index[i])] = solar_wind.calculate_MW(cutout, wind_assets.loc[wind_assets.index[i]], 'wind')

    for _,row in solar_assets.iterrows():
        if row['asset_id'] not in wind_gen_dict.keys(): 
            wind_gen_dict[row['asset_id']] = solar_wind.calculate_MW(cutout, row, 'solar').squeeze()
        else:
            wind_gen_dict[row['asset_id']] += solar_wind.calculate_MW(cutout, row, 'solar').squeeze()

    pv_generation = pd.DataFrame(wind_gen_dict)


    return pv_generation

#Does some input verification before generating the time series
def main():
    # Load configuration files
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)


    #Try reading the arguments passed in the terminal
    assets_path = cfg['solar']['asset_path']
    cutout_path = utils.get_cutout_path(cfg) # Path(sys.argv[2])
    calibration_flag = cfg['solar']['calibration'] # 0 for no calibration, 1 for calibration
    output_path = cfg['solar']['ts_path']


    #Load the solar_assets
    assets = pd.read_csv(assets_path)

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
    