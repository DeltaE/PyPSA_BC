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
from typing import Optional
from pathlib import Path
import requests

def print_update(level: int=None,
                 message: str="--",
                 alert:Optional[bool]=False):
    """Backward-compatible wrapper that routes to the staged logger.

    level 1 -> INFO (file), level 2/3 -> DEBUG (file), level 100 -> ERROR,
    alert -> WARNING. Console shows only WARNING+; full detail is in the log.
    """
    from pypsa_bc.reporting.logger import get_logger
    log = get_logger()
    if alert:
        log.warning(message)
    elif level == 100:
        log.error(message)
    elif level == 1:
        log.info(message)
    else:
        log.debug(message)
    
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
    

def get_region_polygon(geometry):
    '''
    Create an axis-aligned bounding-box polygon around the input geometry.

    The returned polygon is a rectangular extent (minx, miny, maxx, maxy),
    not the original region shape.
    '''
    if len(geometry) == 0:
        raise ValueError('get_region_polygon received empty geometry.')

    west_lon, south_lat, east_lon, north_lat = geometry.total_bounds
    return shapely.geometry.box(west_lon, south_lat, east_lon, north_lat, ccw=True)

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

    Path(file).parent.mkdir(parents=True, exist_ok=True)
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

    cutout.prepare(monthly_requests=True, 
                    concurrent_requests=False,
                    show_progress=True,
                    data_format="netcdf")

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
    skipped = []
    # Assume generators are connected to lowest voltage bus at their given node_code
    for bus in buses:
        # Skip nameless buses (NaN in buses.csv, e.g. unresolved intertie nodes)
        if not isinstance(bus, str) or "_" not in bus:
            skipped.append(bus)
            continue
        node = "_".join(bus.split('_')[1:])
        voltage = int(bus.split('_')[0])
        if node not in bus_dict.keys():
            bus_dict[node] = voltage
        else:
            bus_dict[node] = min(voltage, bus_dict[node])

    if skipped:
        from pypsa_bc.reporting.logger import get_logger
        get_logger("utils").warning(
            f"create_standard_gen_bus_map: skipped {len(skipped)} bus(es) with "
            f"no valid name: {skipped}")

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
    columns_to_add = ['variable_om_costs',  # CAD/MWh # checked with CODERS; EL_20260724
                      'average_fuel_price_CAD_per_MMBtu', # checked with CODERS; EL_20260724
                    #   'average_fuel_price_CAD_per_GJ', # checked with CODERS; EL_20260724  # Deprecated EL_20260724due to CODERS update
                      'carbon_emissions', # _tCO2eq_per_MWh # checked with CODERS; EL_20260724
                      'heat_rate', # _MMBtu_per_MWh # checked with CODERS; EL_20260724
                      'spinning_reserve_capability',
                      'economic_life', # checked with CODERS; EL_20260724
                      'capital_cost_CAD_per_kW',# checked with CODERS; EL_20260724
                    #   'capital_overhead_CAD_per_kW', # Deprecated EL_20260724due to CODERS update
                    #   'overnight_capital_cost_CAD_per_kW',  # Deprecated EL_20260724 due to CODERS update
                    #   'interest_during_construction_CAD_per_kW', # Deprecated EL_20260724 due to CODERS update
                      'annualized_capital_cost_CAD_per_kwyear', # checked with CODERS; EL_20260724
                      'total_project_cost_CAD_per_kW',# checked with CODERS; EL_20260724
                      ]
    # Callers pass either a CODERS gen_type (e.g. 'solar_PV', 'wind_onshore') or a
    # gen_type_copper code (e.g. 'wind_ons'); match on whichever column holds it.
    mask = gen_generic['gen_type'] == gen_type
    if 'gen_type_copper' in gen_generic.columns:
        mask = mask | (gen_generic['gen_type_copper'] == gen_type)
    if not mask.any():
        raise KeyError(
            f"add_generic_columns: no generation_generic row for '{gen_type}'. "
            f"gen_type values={sorted(gen_generic['gen_type'].dropna().unique())}; "
            f"gen_type_copper values={sorted(gen_generic.get('gen_type_copper', pd.Series(dtype=str)).dropna().unique())}."
        )

    param_dict = gen_generic.loc[mask, columns_to_add].iloc[0].to_dict()

    enriched_wind_assets = assets.assign(**param_dict)

    return enriched_wind_assets

def add_generic_columns_tpp(assets, gen_generic):
    '''
    This functions adds generic columns of costs and efficiencies.
    This is needed before creating the .csv for OSeMOSYS via Otoole.
    '''
    # Current generation_generic column names (Oct 2025 CODERS schema).
    columns_to_add = ['variable_om_costs',              # was variable_om_cost_CAD_per_MWh
                      'average_fuel_price_CAD_per_MMBtu',
                      'carbon_emissions',
                      'heat_rate',
                      'spinning_reserve_capability',
                      'economic_life',
                      'capital_cost_CAD_per_kW',
                      'annualized_capital_cost_CAD_per_kwyear',  # was ..._MWyear
                      'total_project_cost_CAD_per_kW',           # was total_project_cost_2020_...
                      'ramp_rate_percent_per_min',
                      "min_up_time_hours",
                      "min_down_time_hours",
                      ]
    # Defensive: only pull columns that exist, warn on the rest (don't hard-crash
    # on a future CODERS generation_generic change).
    have = [c for c in columns_to_add if c in gen_generic.columns]
    missing = [c for c in columns_to_add if c not in gen_generic.columns]
    if missing:
        from pypsa_bc.reporting.logger import get_logger
        get_logger("utils").warning(
            f"add_generic_columns_tpp: generation_generic missing {missing}; skipped")

    df_new = pd.merge(assets, gen_generic[have + ['gen_type']], on='gen_type', how='left')

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