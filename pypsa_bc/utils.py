import yaml
import os
import atlite
import shapely
import geopandas as gpd
import pypsa
import numpy as np
import pandas as pd
from shapely.ops import unary_union
import sys
from pathlib import Path
from colorama import Fore, Style
from typing import Optional, List, Dict
from pathlib import Path
import requests

def print_update(level: int=None,
                 message: str="--",
                 alert:Optional[bool]=False):
    if level is not None:
        if level == 1:
            color = Fore.YELLOW
            prefix="└"
        elif level == 2:
            color = Fore.CYAN
            prefix=" └"
        elif level > 2 and level<100:
            color = Fore.LIGHTBLACK_EX + Style.DIM
            prefix="  └─"
        elif level==100:
            color = Fore.RED
            prefix="└"
        elif alert:
            level=2
            color = Fore.RED
            prefix=" └ X "
    else:
        color = Fore.LIGHTMAGENTA_EX + Style.DIM
        prefix=" ─"
    
    print(f"{color}{prefix}> {message}{Style.RESET_ALL}")

def load_network(network_path):
    network = pypsa.Network(override_component_attrs=get_multi_link_override())
    pypsa.Network.import_from_netcdf(network=network, path=network_path)
    return network

def get_networks(pypsa_results_path:str|Path="results/pypsa")->list:
    """
    Load all the networks in the folder and return a list of the network names
    
    Args:
        pypsa_results_path (str|Path): The folder containing the networks. Default is set to "results/pypsa"
    Returns:
        network_names (list): A list of the network names
    """
    pypsa_results_path=Path(pypsa_results_path)
    network_names = []
    network_dict = {}

    # Loop over each file in the folder
    for file_name in os.listdir(pypsa_results_path):
        if file_name.endswith('.nc'):  # Process only .nc files      
            # pypsa_n_XXXX_coordinated_scale_YYYYMMDD.nc
            network_name = file_name.rstrip('.nc')
            network_name.split('_')[2]
            network_name.split('_')[3]
            network_name.split('_')[4]
            network_name.split('_')[5]
            
            # Load the network and store it in a dictionary
            network_dict[network_name] = load_network(pypsa_results_path / file_name)
           
            # Apply the load_network function and assign the result to the dynamically created variable
            # globals()[network_name] = load_network(pypsa_results_path/file_name)
            
            network_names.append(network_name)
            # Print the assigned variable name for verification
            print_update(level=3,message=f"Assigned: {network_name} = {pypsa_results_path/file_name}")
    return network_names, network_dict

def ev_load_only(network:pypsa.Network,
                 exclude_load_bus:Optional[str]='BC ELC Load'):
    
    load_data = network.loads_t.p_set
    load_data.index = pd.to_datetime(load_data.index)
    exclude_load_bus = 'BC ELC Load'
    load_data_without_excluded = load_data.drop(columns=[exclude_load_bus], errors='ignore')
    ev_load_only = load_data_without_excluded.sum(axis=1)
    return ev_load_only


def get_ev_reults(network_names:list,
                  network_dict:Dict[str,pypsa.Network])->dict:
    ev_load_results:dict = {}

    # Iterate over all networks in network_names
    for network_name in network_names:
        # Retrieve the network object using globals()
        network = network_dict[network_name]
        
        # Apply the ev_load_only function to the network
        result = ev_load_only(network)
        
        # Store the result with the network name in the dictionary
        ev_load_results[network_name] = result
    
    print_update(level=2,message=  "----Summary of EV loads from Results---")
    # Iterate over the ev_load_results dictionary
    for network_name, ev_load in ev_load_results.items():
        # Calculate the sum of EV load for the current network
        total_ev_load = ev_load.sum()
        
        # Print the network name and the total EV load
        print_update(level=3,message=f"Network: {network_name}, Total EV Load: {int(total_ev_load/1E3)} GWh")
    # Now `ev_load_results` contains the ev_load_only result for each network, keyed by network name
    return ev_load_results

def merge_assets(df,subset,sum_list):
    '''
    Function used to reduce hydroelectric datasets from turbines to an aggregate asset.
    This aggregation operator is currently applied only to the installed capacities for the units and
    the annual_avg_energy.
    df: Dataframe which is passed via a groupby operation.
    subset: Name of columns to use for deduplication
    sum_list: Name of parameters/columns to aggegtate using the sum operation.
    Example:
    Input dataframe has following entries below:
    component_id | asset_id | capacity | annual_avg_energy
    BC_MCA01_GEN | BC_MCA_GSS | 492 | 1936.79
    BC_MCA02_GEN | BC_MCA_GSS | 492 | 1936.79
    BC_MCA03_GEN | BC_MCA_GSS | 494 | 1942.7
    BC_MCA04_GEN | BC_MCA_GSS | 494 | 1942.7
    BC_MCA05_GEN | BC_MCA_GSS | 500 | 1968.45
    BC_MCA06_GEN | BC_MCA_GSS | 500 | 1968.45

    Output dataframe will have the following:
    asset_id | capacity | annual_avg_energy
    BC_MCA_GSS | 2972 | 11695.88

    '''
    # Other columns don't matter here for calcualting inflow and associated power production.
    df_out = df.drop_duplicates(subset=subset).set_index("connecting_node_code").copy()
    for param in sum_list:
        df_out[param] = df[param].sum()
    return df_out
        
def set_root(current_dir = Path.cwd()):
    # Check if the last folder name matches 'BC_Combined_Modelling'
    if current_dir.name == 'BC_Combined_Modelling':
        print("The current directory is 'BC_Combined_Modelling'")
    else:
        # Set it as the working directory
        os.chdir('..')

        # Verify the change
        print(f"Current working directory set to: {os.getcwd()}")  

def find_project_root(marker="bc_combined_modelling", start_dir: str = None) -> Path:
        """
        Traverse upwards to find the root directory containing a specific folder or file.

        Args:
            marker (str): The folder or file name to identify the project root.
            start_dir (str): The starting directory for the search. Defaults to current directory.

        Returns:
            Path: The path to the project root.

        Raises:
            FileNotFoundError: If the marker is not found.
        """
        start_dir = Path(start_dir or os.getcwd())
        current_dir = start_dir.resolve()
        while current_dir != current_dir.parent:  # Stop at the filesystem root
            if (current_dir / marker).exists():  # Check for the marker
                print(f"Project root found at: {current_dir}")
                return current_dir
            current_dir = current_dir.parent
        
        raise FileNotFoundError(f"Project root with marker '{marker}' not found starting from '{start_dir}'.")


def setup_environment(marker="bc_combined_modelling", start_dir: str = None):
        """
        Configure the environment for running scripts or notebooks from any directory.

        Args:
            marker (str): The folder or file name to identify the project root.
            start_dir (str): The starting directory for the search. Defaults to current directory.
        """
        try:
            project_root = find_project_root(marker, start_dir=start_dir)
            if str(project_root) not in sys.path:
                sys.path.append(str(project_root))  # Add project root to sys.path
            print(f"Environment set up. Project root added to sys.path: {project_root}")
        except FileNotFoundError as e:
            print(str(e))
            raise

def download_data(source_URL: str, file_path: str) -> str:
    """
    Downloads a file from a given URL and saves it to the specified file path.
    
    Parameters:
        source_URL (str): URL of the file to download.
        file_path (str): Path where the downloaded file will be saved.
    
    Returns:
        str: The file path if download is successful; otherwise, an instruction message.
    """
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    print_update(level=2,message=f"Downloading datafile from: {source_URL}")
    try:
        # Send HTTP GET request
        response = requests.get(source_URL, headers=headers, timeout=30)
        
        # Check if the request was successful
        if response.status_code == 200:
            with open(file_path, 'wb') as file:
                file.write(response.content)
            print_update(level=3,message=f"File downloaded successfully and saved as {file_path}")
            return file_path
        else:
            print_update(level=2,message=f"Failed to download the file. Status code: {response.status_code}")
            return print_update(level=1,message=f"Please download the data manually from {source_URL} and save it to {file_path}")
    except requests.RequestException as e:
        print_update(level=2,message=f">> An error occurred while downloading the file: {e}")
        return f">> Please download the data manually from {source_URL} and save it to {file_path}"

# In a Jupyter Notebook
# Add the following at the start of your notebook:
# 
# from your_script_name import ProjectUtils
# ProjectUtils.setup_environment()

def parse_data_value(value):
    """
    Attempts to convert strings that represent lists (e.g. '[100, 200, 300]')
    into actual lists. If not possible, returns the value as-is.
    """
    if isinstance(value, str):
        try:
            value = eval(value)
        except (SyntaxError, NameError):
            # If eval fails, keep the original string
            pass
    return value


def load_config(config_file):
    '''
    This function loads the configuration file for PyPSA_BC
    config_file: Path + filename of the configuration file. (i.e. ../config/config.yaml)
    '''
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)
    return cfg

def check_path(dir=str|Path):
    '''
    This wrapper handles missing folders and creates them if they do not exist.
    '''
    Path(dir).mkdir(parents=True, exist_ok=True)
    

def get_region_polygon(geometry):
    '''
    This function finds a bounding box of the region and creates a polygon for it.
    Returns a polygon of the regions max/min bounds in terms of lats and lons.
    '''
    if len(geometry) == 1:
        west_lon= geometry.bounds['minx'].iloc[0]
        south_lat  = geometry.bounds['miny'].iloc[0]
        east_lon = geometry.bounds['maxx'].iloc[0]
        north_lat = geometry.bounds['maxy'].iloc[0]
        bbox = (west_lon, south_lat, east_lon, north_lat)
        polygon = shapely.geometry.box(*bbox, ccw=True)
    else:
        print('There remains multiple geometries')
        exit(1)
    
    return polygon

def get_bounds(polygon_list):
    '''
    This function takes in a list of polygons and returns the maximum bounds for them.
    '''
    bounds = {}
    merged_poly = gpd.GeoSeries(unary_union(polygon_list))
    bounds["west_lon"] = merged_poly.geometry.bounds['minx'].iloc[0]
    bounds["south_lat"] = merged_poly.geometry.bounds['miny'].iloc[0]
    bounds["east_lon"] = merged_poly.geometry.bounds['maxx'].iloc[0]
    bounds["north_lat"] = merged_poly.geometry.bounds['maxy'].iloc[0]
    
    return bounds

def get_cutout_path(cfg):
    '''
    This function return the unique name based on the region and start/end year
    for a cutout. 
    return: file path + name for the cutout described by selections in the
    cutout configuration.
    '''
    # Create file_path name with custom year_+date
    start_year = cfg['cutout']["snapshots"]["start"][0][:4]
    end_year = cfg['cutout']["snapshots"]["end"][0][:4]
    prefix = cfg['cutout']['root'] + "BC"  # automate the region name
    if start_year == end_year:
        suffix = start_year
        file = "_".join([prefix, suffix + ".nc"])
    else: # multi_year_file
        suffix = "_".join([start_year, end_year])
        file = "_".join([prefix, suffix + ".nc"])

    return file

def create_era5_cutout(bounds, cfg):
    '''
    This function creates a cutout based on data for era5.
    '''
    # Extract parameters from configuration file
    dx,dy = cfg["cutout"]["dx"], cfg["cutout"]['dy']
    time_horizon = slice(cfg["cutout"]["snapshots"]['start'][0],
                        cfg["cutout"]["snapshots"]['end'][0])

    # get path + filename for the cutout
    file = get_cutout_path(cfg)

    # Create the cutout based on bounds found from above
    cutout = atlite.Cutout(path=file,
                    module=cfg["cutout"]["module"],
                    x=slice(bounds['west_lon'] - dx, bounds['east_lon'] + dx),
                    y=slice(bounds['south_lat'] - dy, bounds['north_lat'] + dy ),
                    dx=dx,
                    dy=dy,
                    time=time_horizon)

    cutout.prepare()

def convert_cid_2_aid(cid,old_aid):
    '''
    This creates an asset id (aid) based on the component id (cid).
    Common Example: 
            cid -> BC_ZBL03_GEN
            old_aid -> BC_ZBL_GSS
            new_aid -> BC_ZBL_GSS
           Example:
           cid -> BC_BR0101_GEN
           old_aid -> BC_BR1_DFS
           new_aid -> BC_BR1_DFS
           Example:
           cid -> BC_BR0102_GEN
           old_aid -> BC_BR2_GSS
           new_aid -> BC_BR2_GSS

    '''
    aid_start = old_aid.split('_')[0]
    cid_2_aid = cid.split('_')[1][:3]
    aid_end = old_aid.split('_')[-1]
    
    if aid_start != cid.split('_')[0]: # error check
        print('Error detected in convert_cid_2_aid')
        exit(3)
    new_aid= "_".join([aid_start, cid_2_aid, aid_end])
    return new_aid 

def write_pickle(data_dict, filepath):
    '''
    Write a pickle file based on a dictionary.
    '''
    with open(filepath,"wb") as f:
        pd.to_pickle(data_dict, f)
    f.close()

def read_pickle(filepath):
    '''
    Read a json file based on a dictionary.
    '''
    with open(filepath, 'rb') as f:
        data_dict = pd.read_pickle(f) 
    return data_dict

def create_standard_gen_bus_map(buses):
    '''
    This function accepts a list a buses and returns a mapping from the buses to the lowest voltage
    bus for each unique node code. The underlying assumption used in selecting a bus to connect a generator to is that
    generators are connected the lowest voltage bus at their given node location.
    Example:
    Buses = ["230_ABN_GSS", "138_ABN_GSS","500_MCA_GSS", "63_MCA_GSS"] -> bus_dict = {ABN_GSS:138, MCA_GSS:63}
    '''
    bus_dict = {}
    # Assume generators are connected to lowest voltage bus at their given node_code
    for bus in buses:
        node = "_".join(bus.split('_')[1:])
        voltage = int(bus.split('_')[0])
        if node not in bus_dict.keys():
            bus_dict[node] = voltage
        else:
            bus_dict[node] = min(voltage,bus_dict[node])

    return bus_dict

def get_gen_bus(node_code, bus_dict):
    '''
    This function returns the correct standardized electric bus for a generator.
    Example:
    node_code = "BC_ABS_GSS" 
    bus_dict = {"ABS_GSS":230, "MCA_GSS":63}
    return -> "230_ABS_GSS"
    '''

    node_code_suffix = "_".join(node_code.split('_')[1:]) # i.e.  BC_ABN_GSS -> ABN_GSS
    return "_".join([str(bus_dict[node_code_suffix]), node_code_suffix])

def get_multi_link_override():
    """
    Gets the multi-link override configuration for PyPSA components.

    This function modifies the default component attributes of PyPSA to allow 
    a single link to have two outputs, which is necessary for modeling cascaded 
    hydroelectric systems. It adds attributes for a second bus, its efficiency, 
    and its output power to the 'Link' component.

    Returns:
        dict: A dictionary with the overridden component attributes for PyPSA.
    """

    # From PyPSA CHP Example: This ensures we can add 2 outputs for a single link i.e bus0 -> bus_1 AND bus_2
    # override_component_attrs = pypsa.descriptors.Dict(
    #     {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    # )
    override_component_attrs = pypsa.components.component_attrs.copy() # new version compatibility
    
    override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "2nd bus",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        1.0,
        "2nd bus efficiency",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "2nd bus output",
        "Output",
    ]

    
    print_update(level=2,message='For more about custom components of pypsa, see >> https://pypsa.readthedocs.io/en/latest/user-guide/components.html#custom-components')
    return override_component_attrs

def get_multi_gen_override():
    '''

    '''
    override_component_attrs = pypsa.descriptors.Dict(
        {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    )
    override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "2nd bus",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        1.0,
        "2nd bus efficiency",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "2nd bus output",
        "Output",
    ]


def add_generic_columns(assets, gen_generic, gen_type):
    '''
    This functions adds generic columns of costs and efficiencies.
    This is needed before creating the .csv for OSeMOSYS via Otoole.
    # MODIFIED 2024-10-01 for CODERS update
    '''
    columns_to_add = ['variable_om_cost_CAD_per_MWh',
                      'average_fuel_price_CAD_per_MMBtu',
                      'average_fuel_price_CAD_per_GJ',
                      'carbon_emissions', # _tCO2eq_per_MWh
                      'heat_rate', # _MMBtu_per_MWh
                      'spinning_reserve_capability',
                      'economic_life',
                      'capital_cost_CAD_per_kW',
                      'capital_overhead_CAD_per_kW',
                      'overnight_capital_cost_CAD_per_kW',
                      'interest_during_construction_CAD_per_kW',
                      'annualized_capital_cost_CAD_per_MWyear',
                      'total_project_cost_2020_CAD_per_kW',
                      ]
    mask = gen_generic['gen_type'] == gen_type # "Wind_onshore"

    param_dict = gen_generic.loc[mask,columns_to_add].iloc[0].to_dict()

    enriched_wind_assets = assets.assign(**param_dict)

    return enriched_wind_assets

def add_generic_columns_tpp(assets, gen_generic):
    '''
    This functions adds generic columns of costs and efficiencies.
    This is needed before creating the .csv for OSeMOSYS via Otoole.
    '''
    columns_to_add = ['variable_om_cost_CAD_per_MWh',
                      'average_fuel_price_CAD_per_MMBtu',
                      'average_fuel_price_CAD_per_GJ',
                      'carbon_emissions',
                      'heat_rate',
                      'spinning_reserve_capability',
                      'economic_life',
                      'capital_cost_CAD_per_kW',
                      'capital_overhead_CAD_per_kW',
                      'overnight_capital_cost_CAD_per_kW',
                      'interest_during_construction_CAD_per_kW',
                      'annualized_capital_cost_CAD_per_MWyear',
                      'total_project_cost_2020_CAD_per_kW',
                      'ramp_rate_percent_per_min',
                      "min_up_time_hours",
                      "min_down_time_hours",
                      ]
    
    # Need to find the matching "generation_type" based on the 
    # assets = assets.rename(columns={'gen_type': 'generation_type'}) 
    df_new = pd.merge(assets, gen_generic[columns_to_add + ['gen_type']], on='gen_type', how='left')


    return df_new

def fix_df_ts_index(
    df:pd.DataFrame, 
    start_date:str='2021-01-01 00:00:00', 
    end_date:str='2021-12-31 23:00:00'):
    '''
    This function hardcodes and fixes the timeseries to be an 8760 timeseries
    beginning in 2021-01-01.
    '''
    new_indices = pd.date_range(start = start_date, end = end_date, freq='h')
    
    df.index = new_indices

    return df

def fix_coders_update(data,col_to_correct,codes):
    '''
    Modified 2024-10-01: CODERS updates to buses creating mismatches in string match incorrectly.
    data:
    col_to_correct: column name containing the values in which we plan to map to corrections for CODERS.
    codes: cords to map from an old string to a new one, example of codes...
    codes = {'BC_BR0101_GEN':'BC_BR101_GEN','BC_BR0102_GEN':'BC_BR102_GEN','BC_BR0103_GEN':'BC_BR103_GEN','BC_BR0104_GEN':'BC_BR104_GEN',
                   'BC_BR0201_GEN':'BC_BR201_GEN','BC_BR0202_GEN':'BC_BR202_GEN','BC_BR0203_GEN':'BC_BR203_GEN','BC_BR0204_GEN':'BC_BR204_GEN'}
    '''
    mask = data[col_to_correct].isin(codes)
    data.loc[mask,col_to_correct] = data.loc[mask,col_to_correct].apply(lambda old_cold: codes[old_cold])

def get_network_features(network_names:str,
                         network_dict:Dict[str,pypsa.Network]):
    # Initialize a dictionary to store results for each network
    network_features = {}

    # Iterate over each network variable name in network_names
    for var_name in network_names:
        # n = globals()[var_name]  # Access the loaded network using the variable name
        n = network_dict[var_name] 
        # Identify relevant columns for each type of feature
        res_link_cols = n.links.index[n.links.index.str.contains('Discharge Link')]
        ror_gen_cols = n.generators.index[n.generators.index.str.contains('RoR')]
        wind_gen_cols = n.generators.index[n.generators.index.str.contains('Wind')]
        solar_gen_cols = n.generators.index[n.generators.index.str.contains('Solar')]
        discharge_cols = n.links.index[n.links.index.str.contains('Discharger')]
        charge_cols = n.links.index[n.links.index.str.contains('Charger')]
        new_wind_gen_cols = n.generators.index[n.generators.index.str.contains('New Wind')]
        new_solar_gen_cols = n.generators.index[n.generators.index.str.contains('New PV')]
        new_backstop_col = n.generators.index[n.generators.index.str.contains('Backstop')]

        # Calculate the sums for each feature
        res_gen = n.links_t.p1[res_link_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        ror_gen = n.generators_t.p[ror_gen_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        wind_gen = n.generators_t.p[wind_gen_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        solar_gen = n.generators_t.p[solar_gen_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        discharge_gen = n.links_t.p1[discharge_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        charge_gen = n.links_t.p1[charge_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        new_wind_gen = n.generators_t.p[new_wind_gen_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        new_solar_gen = n.generators_t.p[new_solar_gen_cols].apply(lambda x: abs(x), axis=0).sum(axis=1)
        backstop_gen = n.generators_t.p[new_backstop_col].apply(lambda x: abs(x), axis=0).sum(axis=1)

        new_capacity = n.generators['p_nom_opt'] - n.generators['p_nom']  # Difference in optimized and initial capacities
        total_new_installed_capacity = new_capacity.sum()
        # Calculate total generation
        tot_gen = res_gen + ror_gen + wind_gen + solar_gen

        # Store the results in the dictionary
        network_features[var_name] = {
            'res_gen': res_gen,
            'ror_gen': ror_gen,
            'wind_gen': wind_gen,
            'solar_gen': solar_gen,
            'discharge_gen': discharge_gen,
            'charge_gen': charge_gen,
            'new_wind_gen': new_wind_gen,
            'new_solar_gen': new_solar_gen,
            'backstop_gen': backstop_gen,
            'total_gen': tot_gen,
            'total_new_installed_capacity': total_new_installed_capacity
        }
        
    return network_features