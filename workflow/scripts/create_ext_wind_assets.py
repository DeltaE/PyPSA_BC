import pandas as pd
import json
from pathlib import Path
from pypsa_bc import utils

# handles the config loading centrally
from pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()

'''
====================================================================================================

    THE ACTUAL FUNCTIONS THAT PRODUCE wind_assets.csv ARE BELOW HERE

====================================================================================================
'''


# Canonical CWTD schema; on-disk headers are mapped to these case-insensitively
# so the curated export ('Project name') and the raw FGP file ('Project Name')
# both work.
_CWTD_CANON = [
    "OBJECTID", "Province/Territory", "Project name", "Total project capacity (MW)",
    "Turbine identifier", "Turbine number in project", "Turbine rated capacity (kW)",
    "Rotor diameter (m)", "Hub height (m)", "Manufacturer", "Model",
    "Commissioning date", "Latitude", "Longitude", "Notes",
]
_CWTD_BY_LOWER = {c.lower(): c for c in _CWTD_CANON}


def _canon_col(c) -> str:
    base = str(c).strip().split(" / ")[0].strip()   # drop bilingual ' / Nom...' suffix
    return _CWTD_BY_LOWER.get(base.lower(), base)    # map any casing to canonical


def load_cwtd(path):
    """Load the Canadian Wind Turbine Database robustly (curated or raw FGP).

    Scans sheets for the turbine table and canonicalises headers (whitespace,
    case, bilingual). Raises a clear error naming the columns if not found.
    """
    xl = pd.ExcelFile(path)
    for sheet in xl.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        df = df.rename(columns={c: _canon_col(c) for c in df.columns})
        if {"Project name", "Model", "Manufacturer"}.issubset(df.columns):
            return df
    first = pd.read_excel(path, sheet_name=xl.sheet_names[0])
    raise KeyError(
        f"CWTD at {path}: no sheet with 'Project name'/'Model'/'Manufacturer'. "
        f"Sheets={xl.sheet_names}; first-sheet columns="
        f"{[str(c).strip() for c in first.columns]}."
    )


#This part builds the wind_assets data frame to be written into a CSV file
#Used in main()
#coders = CODERS data frame from generators.csv
#canada_turbines = Canadian wind turbine data frame from the .xlsx file
#turbine_dict = Dictionary for OEDB configs from a manually written .json file
def generate_wind_assets(wind_assets, turbines, turbine_dict, province):
    '''
    This function creates the wind ...

    NOTE: Grouse turbine removed since it is not operated.
    '''
    #Filter to wind generators in BC for now
    # wind_bc = coders.loc[coders.province.eq('BC') & coders.gen_type_copper.eq('wind')]# [['gen_node_code', 'project_name', 'location', 'latitude', 'longitude', 'annual_avg_energy_unit_in_gwh/y']]
    wind_assets = wind_assets.rename(columns={'annual_avg_energy_unit_in_gwh/y': 'CODERS AAG (GWh/y)'})
    

    # NOTE: Regions need a better handling mechanism

    #Rename some projects for the merge later below
    # BC Renaming
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Dokie'].index.values[0], 'project_name'] = 'Dokie Ridge Wind Farm'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Quality'].index.values[0], 'project_name'] = 'Quality Wind'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Shinish Creek'].index.values[0], 'project_name'] = 'Shinish'

    # Remove the Grouse Mountain turbine
    # turbines = turbines.drop(turbines[turbines['Project name'] == 'Grouse Mountain'].index)

    # Alberta Renaming
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Ardenville'].index.values[0], 'project_name'] = 'Ardenville Wind Farm'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Blackspring Ridge'].index.values[0], 'project_name'] = 'Blackspring Ridge Wind'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Castle River'].index.values[0], 'project_name'] = 'Castle River Wind Farm'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Castlerock Ridge I'].index.values[0], 'project_name'] = 'Castle Rock Ridge'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Chin Chute Wind'].index.values[0], 'project_name'] = 'Chin Chute'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Ghost Pine'].index.values[0], 'project_name'] = 'Ghost Pine Wind Farm'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Halkirk'].index.values[0], 'project_name'] = 'Halkirk Wind'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Oldman Wind'].index.values[0], 'project_name'] = 'Oldman River'
    # wind_assets.at[wind_assets.loc[wind_assets['project_name'] == 'Sinnott'].index.values[0], 'project_name'] = 'Sinnott Wind Farm'
    

    # Drop Optimist Wind Energy
    # turbines = turbines.drop(turbines[turbines['Project name'] == 'Optimist Wind Energy'].index)


    # Saskatchewan Renaming

    # turbines = turbines.drop(turbines[turbines['Project name'] == 'Optimist Wind Energy'].index)


    #Create a combined data frame that is nearly in the final state to be written to wind_assets.csv
    # wind_assets = unique_models.merge(wind_assets.rename(columns={'project_name': 'Project name'}))


    wind_assets = wind_assets.rename(columns={'project_name': 'Project name'})
    # Remove the Grouse Mountain turbine
    turbines = turbines.drop(turbines[turbines['Project name'] == 'Grouse Mountain'].index)
    # Drop Optimist Wind Energy
    turbines = turbines.drop(turbines[turbines['Project name'] == 'Optimist Wind Energy'].index)
    # Loop over each wind asset in the CODERS dataset and see if there is a matching turbine in the canada wind dataset
        # Match: Add details for the matching turbine
        # No Match: 
    # Filter down to just the turbines in BC
    new_wind = []
    for _,row in wind_assets.iterrows():
        row["asset_id"] = row["generation_unit_code"]
        row['Install capacity'] = row["facility_installed_capacity"]
        if province == 'AB': # NOTE: Temporary fix
            mask = (turbines['Project name'] == row['generation_unit_name']) & (row['generation_unit_name'] != 'Cypress') 
        else:
            mask = (turbines['Project name'] == row['generation_unit_name'])
        if turbines[mask].shape[0] >= 1 :
            turbine = turbines[mask]
            row['Model'] = turbine_dict[turbine['Model'].mode().iloc[0]]
            row['Manufacturer'] = turbine['Manufacturer'].mode().iloc[0]
            row['Turbine Rated Capacity (kW)'] = turbine['Turbine rated capacity (kW)'].mode().iloc[0]
            row['Rotor Diameter (m)'] = turbine['Rotor diameter (m)'].mode().iloc[0]
            row['Hub Height (m)'] = turbine['Hub height (m)'].mode().iloc[0]
            row['config_oedb'] = turbine_dict[turbine['Model'].mode().iloc[0]]
        else:
            # Default solution
            row['Model'] = "V100/1800"
            row['Manufacturer'] = "Vestas"
            row['Turbine Rated Capacity (kW)'] = 1800
            row['Rotor Diameter (m)'] = 100
            row['Hub Height (m)'] = 100
            row['config_oedb'] = "V100/1800*134"
        
        new_wind.append(row)
    
    assets = pd.DataFrame(new_wind)
    assets.reset_index(drop=True, inplace=True)


    # turbines = canada_turbines.loc[canada_turbines['Province/Territory'].apply(lambda x: x in regions)]


    # #Set up data frame containing all the turbines to search for in the turbine database
    # important_cols = ['Project name', 'Model', 'Manufacturer', 'Turbine rated capacity (kW)',
    #                    'Rotor diameter (m)', 'Hub height (m)'] #Columns to use for search
    # model_names = turbines['Model'].unique() #Models of turbines (most important)


    # #Get the frame of unique turbine models
    # unique_models = turbines[important_cols].loc[turbines['Model'].isin(model_names)].copy().drop_duplicates()

    # #Column of model/id pairs that are readily formatted for OEDB search function
    # unique_models['config_oedb'] = unique_models['Model'].apply(lambda x: turbine_dict[x])


    # #Get turbine counts to calculate install capacity
    # counts = turbines[important_cols].loc[turbines['Model'].isin(model_names)].groupby(['Project name', 'Model']).count().reset_index()
    # counts = counts.rename(columns={'Manufacturer': 'Counts'}).drop(columns=['Turbine rated capacity (kW)', 'Rotor diameter (m)', 'Hub height (m)'])

    # #Merge the counts back into the searchFrame
    # unique_models = unique_models.merge(counts)

    # #Calculate install capacity using (turbine's rated power capacity) * (# of turbines at a location)
    # unique_models['Install capacity'] = unique_models.apply(lambda x:
    #                                         wind.get_power_cap(x['config_oedb']) * x['Counts'], axis=1)



    # #Generating unique component_id and asset_id columns, as well as a 'Flag' column to be used in create_wind_ts.py later
    # #Empty frame to hold the component_id and asset_id columns
    # id_frame = pd.DataFrame(columns=['component_id', 'asset_id'])

    # #Loop through unique generator codes
    # for code in wind_assets['gen_node_code'].unique():
    #     #Get list of indices, [0, 1, 2, ...] for each project
    #     asset_indexer = wind_assets.loc[wind_assets['gen_node_code'] == code].reset_index().index.values + 1

    #     #Empty lists to contain data for to_id_frame
    #     component_id = []
    #     # asset_id = []
    #     # flag = []

    #     #Fill the component_id and asset_id lists
    #     for index in asset_indexer:
    #         component_id.append(code[3:6])
    #         # asset_id.append(code[3:7] + str(index))
    #         # if index == 1:
    #         #     flag.append(1)
    #         # else:
    #         #     flag.append(0)
        
    #     #Data frame to concatenate to id_frame above
    #     to_id_frame = pd.DataFrame(data={'component_id': component_id})

    #     #Concatenate into id_frame
    #     id_frame = pd.concat([id_frame, to_id_frame])

    # #Concatenate the id_frame into wind_assets
    # id_frame = id_frame.reset_index().drop(columns=['index'])
    # wind_assets_final = pd.concat([id_frame, wind_assets], axis=1)

    # # temp fix: asset_id set based on connecting_node_code
    # wind_assets_final["asset_id"] = wind_assets_final["connecting_node_code"]


    #Return the finalized wind_assets data frame to be written into a csv file
    return assets


# Does some input verification before generating the assets
def main():

    cfg = pypsa_aparser.data_cfg                       # data.yaml = paths/dirs only
    regions_cfg = pypsa_aparser.base_network_cfg['regions']   # regions moved to base_network.yaml

    utils.print_update(level=1,message="Preparing existing wind assets...")


    #Try reading the arguments passed in the terminal
    coders_path = cfg["data"]['coders']['generators']
    coders_generic_path = cfg['data']['coders']['gen_generic']
    canada_turbine_path = cfg["data"]['wind']['can_turbines']
    turbine_dict_path = cfg["data"]['wind']['turbine_dict']
    output_path = cfg["output"]['create_ext_wind_assets']['fname']
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    #Try loading in the CSV file and XLSX file into Pandas data frames
    coders = pd.read_csv(coders_path)
    canada_turbine = load_cwtd(canada_turbine_path)   # canonicalises headers (curated or FGP)

    #Try loading in the JSON file as a dictionary
    with open(turbine_dict_path) as f:
        turbine_dict = json.load(f)
        f.close()

    prov_code_2_name =  {'BC': "British Columbia", "AB": "Alberta",
                         'SK':'Saskatchewan', 'MB':'Manitoba'}


    wind_df = coders.loc[(coders['province'].apply(lambda x: x in regions_cfg))
                          & coders.gen_type_copper.eq('wind_ons')]


    regions = set([prov_code_2_name[region] for region in regions_cfg
               if region in prov_code_2_name.keys()])

    #All is good, start generating wind_assets.csv
    province =  regions_cfg[0]
    wind_assets = generate_wind_assets(wind_df, canada_turbine, turbine_dict, province)

    # Add in CODERS data for the wind asssets
    gen_generic = pd.read_csv(coders_generic_path)
    to_wind_assets_csv = utils.add_generic_columns(wind_assets, gen_generic, gen_type="wind_onshore") # gen_type code synced with CWTD, EL_20260724

    #Write wind_assets.csv
    to_wind_assets_csv.to_csv(output_path, index=False)
    utils.print_update(level=2,message=f"Existing wind assets prepared and saved to {output_path}")

    #Return code 0 is for when everything runs without a problem
    return 0

if __name__ == '__main__':
    main()