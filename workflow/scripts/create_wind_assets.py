import sys
from pathlib import Path
import pandas as pd
import json
import atlite

'''
====================================================================================================

    HELPER FUNCTIONS BELOW HERE, TO BE MOVED TO A MODULE EVENTUALLY

====================================================================================================
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
====================================================================================================

    THE ACTUAL FUNCTIONS THAT PRODUCE wind_assets.csv ARE BELOW HERE

====================================================================================================
'''

#This part builds the wind_assets data frame to be written into a CSV file
#Used in main()
#coders = CODERS data frame from generators.csv
#canada_turbines = Canadian wind turbine data frame from the .xlsx file
#turbine_dict = Dictionary for OEDB configs from a manually written .json file
def generate_wind_assets(coders, canada_turbines, turbine_dict):
    #Filter down to just the turbines in BC
    bc_turbines = canada_turbines.loc[canada_turbines['Province/Territory'] == 'British Columbia']

    #Remove the Grouse Mountain turbine
    bc_turbines = bc_turbines.drop(bc_turbines[bc_turbines['Project name'] == 'Grouse Mountain'].index)

    #Set up data frame containing all the turbines to search for in the turbine database
    important_cols = ['Project name', 'Model', 'Manufacturer', 'Turbine rated capacity (kW)', 'Rotor diameter (m)', 'Hub height (m)'] #Columns to use for search
    model_names = bc_turbines['Model'].unique() #Models of turbines (most important)


    #Get the frame of unique turbine models
    unique_models = bc_turbines[important_cols].loc[bc_turbines['Model'].isin(model_names)].copy().drop_duplicates()

    #Column of model/id pairs that are readily formatted for OEDB search function
    unique_models['config_oedb'] = unique_models['Model'].apply(lambda x: turbine_dict[x])


    #Get turbine counts to calculate install capacity
    counts = bc_turbines[important_cols].loc[bc_turbines['Model'].isin(model_names)].groupby(['Project name', 'Model']).count().reset_index()
    counts = counts.rename(columns={'Manufacturer': 'Counts'}).drop(columns=['Turbine rated capacity (kW)', 'Rotor diameter (m)', 'Hub height (m)'])

    #Merge the counts back into the searchFrame
    unique_models = unique_models.merge(counts)

    #Calculate install capacity using (turbine's rated power capacity) * (# of turbines at a location)
    unique_models['Install capacity'] = unique_models.apply(lambda x: get_power_cap(x['config_oedb']) * x['Counts'], axis=1)


    #Filter to wind generators in BC for now
    wind_bc = coders.loc[coders.province.eq('BC') & coders.gen_type_copper.eq('wind')][['gen_node_code', 'project_name', 'location', 'latitude', 'longitude', 'annual_avg_energy_unit_in_gwh/y']]
    wind_bc = wind_bc.rename(columns={'annual_avg_energy_unit_in_gwh/y': 'CODERS AAG (GWh/y)'})

    #Rename some projects for the merge later below
    wind_bc.at[wind_bc.loc[wind_bc['project_name'] == 'Dokie'].index.values[0], 'project_name'] = 'Dokie Ridge Wind Farm'
    wind_bc.at[wind_bc.loc[wind_bc['project_name'] == 'Quality'].index.values[0], 'project_name'] = 'Quality Wind'
    wind_bc.at[wind_bc.loc[wind_bc['project_name'] == 'Shinish Creek'].index.values[0], 'project_name'] = 'Shinish'


    #Create a combined data frame that is nearly in the final state to be written to wind_assets.csv
    wind_assets = unique_models.merge(wind_bc.rename(columns={'project_name': 'Project name'}))


    #Generating unique component_id and asset_id columns, as well as a 'Flag' column to be used in create_wind_ts.py later
    #Empty frame to hold the component_id and asset_id columns
    id_frame = pd.DataFrame(columns=['component_id', 'asset_id', 'Flag'])

    #Loop through unique generator codes
    for code in wind_assets['gen_node_code'].unique():
        #Get list of indices, [0, 1, 2, ...] for each project
        asset_indexer = wind_assets.loc[wind_assets['gen_node_code'] == code].reset_index().index.values + 1

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
        to_id_frame = pd.DataFrame(data={'component_id': component_id, 'asset_id': asset_id, 'Flag': flag})

        #Concatenate into id_frame
        id_frame = pd.concat([id_frame, to_id_frame])

    #Concatenate the id_frame into wind_assets
    id_frame = id_frame.reset_index().drop(columns=['index'])
    wind_assets = pd.concat([id_frame, wind_assets], axis=1)


    #Return the finalized wind_assets data frame to be written into a csv file
    return wind_assets

#Does some input verification before generating the assets
def main():
    try:
        #Try reading the arguments passed in the terminal
        coders_path = Path(sys.argv[1])
        canada_turbine_path = Path(sys.argv[2])
        turbine_dict_path = Path(sys.argv[3])
        output_path = Path(sys.argv[4])
    except:
        #Less than 4 arguments, return error code 1
        print('There are inputs missing')
        return 1
    else:
        #Correct number of arguments
        try:
            #Try loading in the CSV file and XLSX file into Pandas data frames
            coders = pd.read_csv(coders_path)
            canada_turbine = pd.read_excel(canada_turbine_path)

            #Try loading in the JSON file as a dictionary
            with open(turbine_dict_path) as f:
                turbine_dict = json.load(f)
                f.close()
        except:
            #An input may be spelled incorrectly, return error code 2
            print('One or more inputs are in the wrong format')
            return 2
        else:
            #All is good, start generating wind_assets.csv
            to_wind_assets_csv = generate_wind_assets(coders, canada_turbine, turbine_dict)

            #Write wind_assets.csv
            to_wind_assets_csv.to_csv(output_path, index=False)

            #Return code 0 is for when everything runs without a problem
            return 0

if __name__ == '__main__':
    main()
    