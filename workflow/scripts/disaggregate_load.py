import pandas as pd
from pypsa_bc import utils
from pathlib import Path
import geopandas as gpd
import pandas as pd
import folium

import plotly.express as px
import numpy as np
from pypsa_bc import utils
from typing import Optional
import warnings

# handles the config loading centrally
from pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()
plot_save_to_root:Path=pypsa_aparser.get_visual_root
    
# Suppress specific warnings
warnings.filterwarnings("ignore", category=FutureWarning)
# The latest version of data doesn't have the "ORG_TYPE" fields to isolate Regional Districts. 
# This names are collected from previous version of data.
CEEI_regional_districts_name = [
    'Alberni-Clayoquot',
    'Bulkley-Nechako',
    'Capital',
    'Cariboo',
    'Central Coast',
    'Central Kootenay',
    'Central Okanagan',
    'Columbia-Shuswap',
    'Comox Valley',
    'Cowichan Valley',
    'East Kootenay',
    'Fraser Valley',
    'Fraser-Fort George',
    'Metro-Vancouver',
    'Kitimat-Stikine',
    'Mount Waddington',
    'Nanaimo',
    'North Okanagan',
    'Northern Rockies',
    'Okanagan-Similkameen',
    'Peace River',
    'Powell River',
    'Skeena-Queen Charlotte',
    'Squamish-Lillooet',
    'Stikine',
    'Strathcona',
    'Sunshine Coast',
    'Thompson-Nicola',
    'Kootenay Boundary',
    'Stikine Region'
]


#Format the hourly load data to account for inconsistencies with how BC Hydro handles daylight savings
def fix_hourly_load(load_bch_raw:pd.DataFrame,
                    year:int):
    
    #Filter down to just the loads
    load_bch_data = pd.to_numeric(load_bch_raw[load_bch_raw.columns[-1]], errors='coerce')
    #Rightmost column in BC Hydro hourly load spreadsheet, to_numeric converts the labels into NaN
    load_bch_kWh = load_bch_data.reset_index().drop(columns=['index'])
    
    load_bch_kWh.columns = ['LOAD']

    #Remove NaN values (the column labels and sometimes the value for DST hour) and 0 values (sometimes the DST hour), then convert from MWh to kWh
    load_bch_MWh = load_bch_kWh.loc[load_bch_kWh['LOAD'] > 0] * 1000

    #Index hourly loads by hourly timestamp
    load_bch_MWh['TIME'] = pd.date_range(start=str(year)+'-01-01 00:00:00', end=str(year)+'-12-31 23:00:00', freq='h')
    load_bch_MWh = load_bch_MWh.set_index('TIME')

    load_bch_MWh=add_normalize_load_data(load_bch_MWh)
    load_bch_MWh.to_csv(f'data/processed_data/load/Hourly_profile_{year}.csv')
    
    plot_hourly_profile(load_bch_MWh)
    
    #Return the hourly load for the entire province for a given year
    return load_bch_MWh

def add_normalize_load_data(load:pd.DataFrame):
   
    load['LOAD_profile_norm_total2hr'] = load['LOAD'] / load['LOAD'].sum()
    load['LOAD_profile_norm_peak2hr'] = load['LOAD'] / load['LOAD'].max()
    
    return load

def visualize_ratios_in_map(ratios:pd.DataFrame,
                            BC_boundary:gpd.GeoDataFrame):

    # Merge boundary GeoDataFrame with ratios on 'REGION' (ratios) and 'Region' (boundary)
    ratios_ = ratios.merge(BC_boundary[['Region', 'geometry']], left_index=True, right_on='Region', how='left')
    ratios_gdf=gpd.GeoDataFrame(ratios_,geometry='geometry')

    # Assuming 'ratios_gdf' is your GeoDataFrame

    # Initialize the base map
    m = folium.Map(location=[53.7267, -127.6476], zoom_start=5)

    # Add the first layer (PROPORTION_RES)
    folium.Choropleth(
        geo_data=ratios_gdf,
        data=ratios_gdf,
        columns=['Region', 'PROPORTION_RES'],
        key_on='feature.properties.Region',
        fill_color='Reds',
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name='Proportion of Residential Load',
        highlight=True,
        name='Residential Load'
    ).add_to(m)

    # Add the second layer (PROPORTION_CSMI)
    folium.Choropleth(
        geo_data=ratios_gdf,
        data=ratios_gdf,
        columns=['Region', 'PROPORTION_CSMI'],
        key_on='feature.properties.Region',
        fill_color='Blues',
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name='Proportion of Commercial, Small and Medium Industries Load',
        highlight=True,
        name='Commercial Load'
    ).add_to(m)

    # Add layer control to toggle between layers
    folium.LayerControl().add_to(m)
    proportions_vis_save_to= plot_save_to_root/'PROVICIAL_2_Regional_loads_ratios.html'
    # Display the map
    m.save(proportions_vis_save_to)

    utils.print_update(level=2,message=f"Provincial load to Regional Districts load ratios' visual saved to : {proportions_vis_save_to}")

#This part disaggregates the hourly load data from BC Hydro into 27 subdivisions of BC

def get_proportions(ceei:pd.DataFrame,
                    proportions_data_path:Path,
                    force_replace:bool=False):
    
    if not proportions_data_path.exists() or force_replace:
        utils.print_update(level=3,message="Preparing Provincial to Regional Proportions data")
        #Get all the regions of BC to loop through
        regions = ceei['ORG_NAME'].unique()

        #Total load for all of BC (Res + CSMI)
        ceei_total = ceei['CONSUMPTION_TOTAL'].sum()

        #Get proportion of annual loads per region
        proportions = pd.DataFrame(columns=['REGION', 'PROPORTION_RES', 'PROPORTION_CSMI'])

        #Create proportions dataframe, do in sorted alphabetical order so that it's easier to look through manually if needed
        for r in sorted(regions.tolist()):
            #Obtain proportion for residential and industrial separately
            proportion_res = ceei.loc[(ceei['ORG_NAME'] == r) & (ceei['SUB_SECTOR'] == 'Res')]['CONSUMPTION_TOTAL'].sum() / ceei_total
            proportion_csmi = ceei.loc[(ceei['ORG_NAME'] == r) & (ceei['SUB_SECTOR'] == 'CSMI')]['CONSUMPTION_TOTAL'].sum() / ceei_total

            #Append to final proportions dataframe
            to_proportions = pd.DataFrame(data={'REGION': [r], 'PROPORTION_RES': [proportion_res], 'PROPORTION_CSMI': [proportion_csmi]})
            to_proportions = to_proportions.dropna(axis=1, how="all").dropna(axis=0, how="all")
            proportions = pd.concat([proportions, to_proportions])

        proportions = proportions.set_index('REGION')


        #Comox Valley & Strathcona combine into Comox-Strathcona, Metro-Vancouver becomes GreaterVancouver (this is to line up with GADM naming convention)
        proportions.loc['Comox Valley'] += proportions.loc['Strathcona']
        proportions = proportions.drop('Strathcona')
        proportions = proportions.rename(index={'Comox Valley': 'Comox-Strathcona', 'Metro-Vancouver': 'Greater Vancouver'})
    
        proportions_data_path.parent.mkdir(parents=True,exist_ok=True)
        proportions.to_csv(proportions_data_path)
        utils.print_update(level=2,message=f"Provincial load to Regional Districts load ratios' saved to :{proportions_data_path}")
    else:
        utils.print_update(level=3,message="Proportions data already exists. Use force_replace=True to overwrite.")
        proportions=pd.read_csv(proportions_data_path)
        proportions=proportions.set_index('REGION')
        
    return proportions
        
def preprocess_ceei_data(ceei_Buildings_eng_file_path:str|Path,
                         profile_data_year:int)->pd.DataFrame:
    """
    Preprocess the Community Energy and Emissions Inventory (CEEI) data for a given year.
    
    # Job:
    This function reads the CEEI data from the source Excel file, processes it to include only the relevant
    information i.e. electricity usage for Residential (RES), Commercial, Small adn Medium Industries (CSMI) sector for all 28 Regional districts, 
    and saves the processed data to a CSV file.
    
    # Args:
        ceei_Buildings_eng_file_path (str | Path): The file path to the CEEI Buildings Excel file.
        profile_data_year (int): The year for which the profile data is being processed.
        
    # Returns:
        pd.DataFrame: A DataFrame containing the processed CEEI data for the specified year.
        Data: Annual electricity usage for RES, CSMI sector, for all 28 Regional districts.
        
    """
        
    if not ceei_Buildings_eng_file_path.exists():
        utils.print_update(level=3,message="CEEI data not found. Downloading from source...")
        # RES, CSMI data
        utils.download_data(source_URL='https://www2.gov.bc.ca/assets/gov/environment/climate-change/data/ceei/bc_utilities_energy_and_emissions_data_at_the_community_level.xlsx',
                    file_path=ceei_Buildings_eng_file_path)

        #TRA data
        # utils.download_data(source_URL='https://www2.gov.bc.ca/assets/gov/environment/climate-change/data/ceei/bc_on_road_transportation_data_at_the_community_level.xlsx',
        #                     file_path=ceei_tra_path)
        
    ceei = pd.read_excel(ceei_Buildings_eng_file_path, sheet_name='Combined') # The source file is an excel file
    utils.print_update(level=2,message=f"Provincial Energy data from Community Energy Inventory (CEEI) loaded from: {ceei_Buildings_eng_file_path}")
    
    if 'ORG_TYPE' not in ceei.columns: # To handle 2022 datafile's missing ORG_TYPE columns
        utils.print_update(level=3,message="ORG_TYPE column not found in the CEEI data. Creating ORG_TYPE column from ORG_NAME")
        # Filter CEEI data to include only rows where ORG_NAME contains 'District'
        ceei = ceei[ceei['ORG_NAME'].isin(CEEI_regional_districts_name)]

        # Add a new column ORG_TYPE with value 'Regional District'
        ceei['ORG_TYPE'] = 'Regional District'

    # CEEI2022 data goes up to 2022, use 2021 if looking at hourly loads 2021 or later
    # We only need data for electricity usage from regional districts for this disaggregation step
    data_year_selection=min(2021, profile_data_year)
    ceei_data = ceei.loc[(ceei.YEAR == data_year_selection) & (ceei.ENERGY_TYPE == 'ELEC') & (ceei.ORG_TYPE == 'Regional District')]
    ceei_data.to_csv(f'data/processed_data/load/CEEI_{data_year_selection}_RD_ELEC.csv')
    
    return ceei_data


def plot_hourly_profile(data:pd.DataFrame,
                        show:Optional[bool]=False):
    utils.print_update(level=2,message="Plotting the normalized load profile for the year...")

    title = f'Normalized Load Profile {data.index[0].year}'
    peak_load_time = data['LOAD'].idxmax()
    peak_load_value = data['LOAD_profile_norm_peak2hr'].max()

    # Function to resample data
    def resample_data(freq):
        return data.resample(freq).mean()

    # Create the initial figure
    fig = px.line(data, x=data.index, y='LOAD_profile_norm_peak2hr', title=title)
    fig.add_scatter(
        x=[peak_load_time],
        y=[peak_load_value],
        mode='lines',
        marker=dict(color='red', size=8),
        name='Peak Load'
    )

    # Customize layout for minimal and clean look
    fig.update_layout(
        template='plotly_white',
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18}
        },
        xaxis_title=None,
        yaxis_title="Normalized Load Profile",
        xaxis=dict(showgrid=False, ticks="outside"),
        yaxis=dict(showgrid=True, gridcolor="lightgrey"),
        font=dict(size=14),
        margin=dict(l=40, r=40, t=50, b=40),
        showlegend=False
    )

    # Add dropdown for resampling
    freq_options = ['h','D', 'W', 'ME']  # Daily, Weekly, Monthly
    updatemenus = [
        {
            "buttons": [
                {
                    "label": f"Sampling: {freq}",
                    "method": "update",
                    "args": [
                        {"x": [resample_data(freq).index], "y": [resample_data(freq)['LOAD_profile_norm_peak2hr']]},
                        {"title": f"Normalized Load Profile (Resampled: {freq})"}
                    ],
                }
                for freq in freq_options
            ],
            "direction": "down",
            "x": 0.8,
            "xanchor": "center",
            "y": 1.15,
            "yanchor": "top",
            "showactive": True,
            "pad": {"r": 10, "t": 10},
        }
    ]

    fig.update_layout(updatemenus=updatemenus)

    # Show and save the figure

    plot_save_to=str(Path('vis') / title) + '.html'
    fig.write_html(plot_save_to)
    utils.print_update(level=3,message=f"plot saved to {plot_save_to}")
    if show:
        fig.show()

    
def disaggregate_load(proportions:pd.date_range,
                      hourly_profile:pd.DataFrame):
    
    #Construct dataframes for hourly load. 28 subdivisions, each with residential and industrial loads
    regional_res = pd.DataFrame()
    regional_csmi = pd.DataFrame()

    # Apply the proportions to the hourly load data to get the disaggregated hourly load for each region in BC
    for region in proportions.index:
        p_res = proportions.loc[region]['PROPORTION_RES']
        p_csmi = proportions.loc[region]['PROPORTION_CSMI']
        
        #Spaces are removed in the column names, now all region names line up with the GADM names
        regional_res[region.replace(' ', '')] = hourly_profile['LOAD'].apply(lambda x: x * p_res)
        regional_csmi[region.replace(' ', '')] = hourly_profile['LOAD'].apply(lambda x: x * p_csmi)

    #Return the dataframes for disaggregated residential loads and disaggregated industrial loads
    return regional_res, regional_csmi

#Does some input verification before generating the regional hourly loads
""" 
USAGE:

python dissagregate_load.py <ARG1> <ARG2> <ARG3> <ARG4>

ARG1 = The CEEI spreadsheet file
ARG2 = The folder which contains the hourly load data from BC Hydro. This folder should contain files named BalancingAuthorityLoad20XX.xls
ARG3 = The year to use for disaggregation
ARG4 = The folder to write the outputs to. This script outputs two CSV files, one for residential load, one for industrial load

EXAMPLE:

'python disaggregate_load.py CEEI_2020.xlsx .\load\ 2015 .\hourly_load_region\'
    
"""


def main(provincial_total_load_MWh:float):
    """
    Main function to disaggregate hourly load data for PyPSA.
    # Parameters:
        load_source (str): The source of load data to use for disaggregation. Default is 'ceei', supports 'nexus' for Nexus BC data.
    # Returns:
    int: Return code 0 if the function runs without any problems.
    
    The function performs the following steps:
        1. Retrieves configuration settings.
        2. Extracts snapshot start and end times.
        3. Loads Community Energy and Emissions Inventory (CEEI) data.
        4. Downloads necessary data files if they do not exist.
        5. Loads and processes hourly load data.
        6. Disaggregates the load data using the specified method.
        7. Saves the disaggregated load data to output files.
    """

    # Get configuration
    cfg=pypsa_aparser.pypsa_cfg
    utils.print_update(level=1,message="Disaggregating hourly load data for PyPSA...")
    
    # >>> (1) Load Snapshot (start_time, end_time)
    (start_time,end_time) = pypsa_aparser.get_snapshot
    utils.print_update(level=2,message=f"Snapshot extracted: {start_time}, {end_time}")
    year =  int(start_time[:4][:4]) # int(sys.argv[3])
    
    # >>> Get the Energy data, Proportions, Hourly Profile datafile paths
    ceei_building_egy_path = Path(cfg['data']["load"]["ceei"]) 
    # ceei_tra_data_path=Path(cfg['data']["load"]["cee_transport"]) # not in use for now
    output_path_res = Path(cfg['output']["disaggregate_load"]["res_path"])
    output_path_csmi = Path(cfg['output']['disaggregate_load']["csmi_path"]) # CSMI - Commercial and Small-medium Industries
    proportions_data_path=Path(cfg['output']['disaggregate_load']["province_2_RD_proportions"])
    hourly_path = cfg['data']["load"]["bch"] +str(year) + ".xls"      # BCH laod profile data
        
    # >>> (2) Pre-process and load the CEEI energy data from sourcefile
    ceei_data:pd.DataFrame=preprocess_ceei_data(ceei_building_egy_path,
                                                year)

    # >>> (3) Get the Proportions of electricity consumption for each region 

    proportions=get_proportions(ceei_data,
                                proportions_data_path,
                                force_replace= False)

    # Hourly load data needs some fixing
    hourly = fix_hourly_load(pd.read_excel(hourly_path), year)
    utils.print_update(level=2,message=f"Hourly load data loaded from: {hourly_path}")

    hourly['LOAD']=provincial_total_load_MWh*hourly['LOAD_profile_norm_total2hr']

    hourly_res, hourly_csmi = disaggregate_load(proportions,
                                                hourly)
    
    BC_boundary = gpd.read_file('data/processed_data/regions/gadm41_Canada_L2_BC.geojson')
    visualize_ratios_in_map(proportions,
                            BC_boundary)
    hourly_res = hourly_res / 1000 # convert from KW-hr to MW-hr
    hourly_csmi = hourly_csmi / 1000 # convert from KW-hr to MW-hr

    # Write to files to the output folder path
    hourly_res.to_csv(output_path_res)
    utils.print_update(level=2,message=f"Hourly residential load saved to: {output_path_res}")
    
    hourly_csmi.to_csv(output_path_csmi)
    utils.print_update(level=2,message=f"Hourly industrial load saved to: {output_path_csmi}")

    #Return code 0 is for when everything runs without a problem
    return 0
    
if __name__ == '__main__':
    main()