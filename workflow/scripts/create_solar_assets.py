from pathlib import Path
import sys
import pandas as pd
from bc_power import utils

#This part builds the solar_assets data frame to be written into a CSV file
#Used in main()
#coders = CODERS data frame from generators.csv
def generate_solar_assets(coders):
    #Filter down to just the solar panels in BC
    solar_assets = coders.loc[coders.province.eq('BC') & coders.gen_type_copper.eq('solar')]
    solar_assets = solar_assets.rename(columns={'install_capacity_in_mw': 'Install capacity','annual_avg_energy_unit_in_gwh/y': 'CODERS AAG (GWh/y)'}).reset_index().drop(columns=['index'])

    #Generating unique component_id and asset_id columns, as well as a 'Flag' column to be used in create_wind_ts.py later
    #Empty frame to hold the component_id and asset_id columns
    id_frame = pd.DataFrame(columns=['component_id', 'asset_id', 'Flag'])

    #Loop through unique generator codes
    for code in solar_assets['gen_node_code'].unique():
        #Get list of indices, [0, 1, 2, ...] for each project
        asset_indexer = solar_assets.loc[solar_assets['gen_node_code'] == code].reset_index().index.values + 1

        #Empty lists to contain data for to_id_frame
        component_id = []
        asset_id = []
        flag = []

        #Fill the component_id and asset_id lists
        for index in asset_indexer:
            component_id.append(code[3:6])
            asset_id.append(code[3:7] + str(index))
            if index == 1:
                flag.append(1)
            else:
                flag.append(0)
        
        #Data frame to concatenate to id_frame above
        to_id_frame = pd.DataFrame(data={'component_id': component_id, 'Flag': flag})

        #Concatenate into id_frame
        id_frame = pd.concat([id_frame, to_id_frame])

    # Concatenate the id_frame into solar_assets
    id_frame = id_frame.reset_index().drop(columns=['index'])
    solar_assets = pd.concat([id_frame, solar_assets], axis=1)

    # temp fix: asset_id set based on connecting_node_code
    solar_assets["asset_id"] = solar_assets["connecting_node_code"]


    #Return the finalized solar_assets data frame to be written into a csv file
    return solar_assets



#Does some input verification before generating the assets
def main():

    # load configuration files
    config_file = r"config/config2.yaml"
    cfg = utils.load_config(config_file)

    #Try reading the arguments passed in the terminal
    coders_path = cfg["data"]['coders']['generators']
    output_path = cfg["output"]['create_solar_assets']['fname']

    # Try loading in the CSV file into a Pandas data frame
    coders = pd.read_csv(coders_path)

    # Start generating solar_assets.csv
    to_solar_assets_csv = generate_solar_assets(coders)

    # Write solar_assets.csv
    to_solar_assets_csv.to_csv(output_path, index=False)

    #Return code 0 is for when everything runs without a problem
    return 0

if __name__ == '__main__':
    main()
    