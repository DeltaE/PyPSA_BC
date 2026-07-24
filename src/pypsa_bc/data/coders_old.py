# for CANADian power system data only.

import os
import sys
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import pandas as pd
import requests
import yaml
from shapely.geometry import Point

PRINT_LEVEL_BASE=3

# # Ensure the script runs from the project root directory
# project_root = Path(__file__).resolve().parent.parent
# if str(project_root) not in sys.path:
#     sys.path.insert(0, str(project_root))
# os.chdir(project_root)

CODERS_CFG_PATH="data/downloaded_data/CODERS/coders_api.yaml"

def load_config(file_path):

    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)

    return data

def load_api_key(file_path=CODERS_CFG_PATH):
    """
    Loads an API key from a configuration file.
    Args:
        file_path (str): The path to the YAML configuration file containing API keys.
                         Defaults to "CODERS/coders_api.yaml".
    Returns:
        str or None: The API key for the default user if specified and available.
                     If no default user is specified or their key is unavailable,
                     returns the first available API key from the configuration.
                     Returns None if no API key is found.
    """
    try:
        api_cfg = load_config(file_path)
        if api_cfg is None:
            print(f"API key file is empty or could not be loaded: {file_path}")
            print("Please create a YAML file at the above path with the following structure:")
            print("""
        api_keys:
          your_username: your_api_key_here
        Default_user: your_username
            """)
            print(f"save the file to : {file_path} and try again.")
            print("Refer to the CODERS API setup guide for more details.")
            return None
    except FileNotFoundError:
            print(f"API key file not found: {file_path}")
            print("Please create a YAML file at the above path with the following structure:")
            print("""
        api_keys:
          your_username: your_api_key_here
        Default_user: your_username
            """)
            print("Refer to the CODERS API setup guide for more details.")
            return None

    default_user = api_cfg.get("Default_user")
    api_keys = api_cfg.get("api_keys", {})

    if default_user:
        api_key = api_keys.get(default_user)
        if api_key:
            return api_key,default_user

    # fallback: try any other API key
    for user, key in api_keys.items():
        if key:
            return key,user

    return None  # or raise an exception



@dataclass
class CODERSData:
    data_config_path: str | Path = None

    """
    Canadian power system data processor using the CODERS API.
    
    This class provides comprehensive access to Canadian power system infrastructure
    data through the CODERS (Canadian Open Data Exchange for Renewable Energy Systems)
    API. It enables retrieval, caching, and processing of transmission lines,
    substations, generators, and other power system components for renewable energy
    integration analysis.
    
    Key Functionality:
    - API-based data retrieval from CODERS database
    - Local data caching and persistence for improved performance
    - Provincial and national data filtering capabilities
    - Geographic data processing with GeoDataFrame support
    - Data validation and error handling for API operations
    
    Data Sources Available:
    - Power generation facilities (generators)
    - Transmission infrastructure (lines, substations)
    - Regional power system characteristics
    - Provincial energy system data
    
    Inherits from:
        AttributesParser: Base class providing configuration management and regional attributes
        
    Attributes:
        coders_data_config (dict): CODERS-specific configuration parameters
        url (str): Base URL for CODERS API endpoints
        api_user (str): API authentication key for CODERS access
        query (str): Formatted query string with authentication
        data_pull (dict): Configuration for data retrieval and storage
        table_list (list): Available data tables from configuration
        region_data (pd.DataFrame/gpd.GeoDataFrame): Filtered regional data
        
    API Requirements:
        - Valid CODERS API key (stored in coders_api.yaml)
        - Network connectivity for data retrieval
        - Proper authentication configuration
        
    Example:
        >>> coders = CODERSData(
        ...     config_file_path="config/config_BC.yaml",
        ...     region_short_code="BC"
        ... )
        >>> 
        >>> # Get provincial transmission data
        >>> bc_substations = coders.get_table_provincial('substations')
        >>> 
        >>> # Get national generator data with forced update
        >>> generators_df, generators_gdf = coders.get_table_canada(
        ...     'generators', 
        ...     force_update=True
        ... )
        
    Data Persistence:
        - Automatic local caching reduces API calls
        - Pickle format for efficient data storage
        - Configurable data refresh policies
        - Regional data filtering and storage
        
    Notes:
        - Requires active internet connection for initial data retrieval
        - API rate limits may apply for excessive requests
        - Local data cache improves performance for repeated analyses
        - Geographic data automatically converted to EPSG:4326 projection
        
    References:
        - CODERS API Documentation: https://sesit.dev/api/docs
        - Canadian power system data standards and formats
    """
    
    def __post_init__(self):
        """
        Initialize inherited attributes and CODERS API configuration.
        
        This method:
        1. Calls parent __post_init__ to inherit configuration and regional attributes
        2. Loads CODERS-specific configuration from config files
        3. Sets up API authentication and connection parameters
        4. Initializes data retrieval and storage configuration
        5. Prepares table list and query formatting
        
        Inherited attributes from AttributesParser:
        - Configuration file parsing and validation
        - Regional identification (region_short_code, region_code_validity)
        - Data storage paths and directory management
        
        CODERS Configuration:
        - API endpoint URLs and authentication
        - Data table specifications and requirements
        - Local storage paths and file naming conventions
        - Regional filtering and validation parameters
        
        Raises:
            ConfigurationError: If CODERS configuration is missing or invalid
            AuthenticationError: If API key is not properly configured
            NetworkError: If API connectivity cannot be established
        """
        

        # Load CODERS data config
        self.config=load_config(self.data_config_path)
        self.api_cfg_path=self.config.get('coders', {}).get('api_config', CODERS_CFG_PATH)
        
        api_key,user = load_api_key(self.api_cfg_path)
        
        if api_key is None:
            print("No API key found. Please ensure you have a valid API key in the configuration file.")
        else:
            print(f"CODERS API key loaded from: {self.api_cfg_path}")
            print(f"user {user}: {api_key}")
        
        
        self.coders_data_config = self.config.get('coders')
        self.url = self.coders_data_config.get('url', None)
        self.api_user=api_key
        
        self.api_query = f"?key={self.api_user}"
        self.data_pull = self.coders_data_config.get('data_pull', {})
        self.table_list = self.coders_data_config.get('tables', [])

    def is_table_name_required(self, table_name: str):
        """
        Validate if a specified table name is configured and required for analysis.
        
        This method checks whether a requested data table is included in the
        configured list of required tables for the current analysis. It serves
        as a validation gate to prevent unnecessary API calls and ensure only
        relevant data is processed.
        
        Args:
            table_name (str): Name of the data table to validate
                            (e.g., 'generators', 'substations', 'transmission_lines')
                            
        Returns:
            bool: True if table is configured and required, False otherwise
            
        Example:
            >>> coders = CODERSData(**config)
            >>> if coders.is_table_name_required('generators'):
            ...     data = coders.get_table_provincial('generators')
        """
        if table_name in self.table_list:
            return True
    
    def show_list(self) -> list:
        """
        Fetch and display available data tables from the CODERS API for a specified source.
        
        This method queries the CODERS API to retrieve and display the complete list
        of available data tables for a given data source. It provides users with
        an inventory of accessible datasets and helps identify appropriate table
        names for data retrieval operations.
        
        Data Sources:
        - 'cef': Canadian Energy Facts data tables
        - 'coders': Core CODERS power system infrastructure tables
        
        Args:
            None
        
        Returns:
            list: List of available table names for the specified source.
                 Returns empty list if API request fails.
        
        Raises:
            RuntimeError: If API returns non-200 status code
            requests.RequestException: If network connectivity issues occur
            
        Example:
            >>> coders = CODERSData(**config)
            >>> 
            >>> # List Canadian Energy Facts tables
            >>> cef_tables = coders.show_list('cef')
            >>> print(f"Available CEF tables: {cef_tables}")
            >>> 
            >>> # List core CODERS infrastructure tables
            >>> coders_tables = coders.show_list('coders')
            >>> print(f"Available CODERS tables: {coders_tables}")
            
        Notes:
            - Requires active internet connection and valid API authentication
            - Table availability may vary based on data source updates
            - Use returned table names for subsequent data retrieval calls
        """
        print(">> Fetching the list of data tables from CODERS")
        try:
            
            response = requests.get(f"{self.url}/tables{self.api_query}")
            if response.status_code == 200:
                tables_list = response.json()
                print(f"CODERS data available:\n {tables_list}")
                return tables_list
            else:
                raise RuntimeError(f">> Error fetching tables list for CODERS: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f">> Connection error while fetching tables list: {e}")
            return []
    
    
    def fetch_data(self, 
                   table_name: str) -> pd.DataFrame:
        """
        Retrieve data from the CODERS API for a specified table.
        
        This method performs direct API calls to fetch power system data from
        the CODERS database. It handles HTTP requests, response validation,
        and data format conversion to return structured pandas DataFrames
        suitable for analysis.
        
        Args:
            table_name (str): Name of the data table to retrieve from CODERS API
                            (e.g., 'generators', 'substations', 'transmission_lines')
                            
        Returns:
            pd.DataFrame: Structured data from the specified CODERS table
            
        Raises:
            RuntimeError: If API returns non-200 status code or request fails
            requests.RequestException: If network connectivity issues occur
            JSONDecodeError: If API response cannot be parsed as valid JSON
            
        Example:
            >>> coders = CODERSData(**config)
            >>> generators_data = coders.fetch_data('generators')
            >>> print(f"Retrieved {len(generators_data)} generator records")
            
        Notes:
            - Requires valid API authentication and network connectivity
            - Raw data retrieval without local caching or persistence
            - Use get_table_canada() or get_table_provincial() for cached access
            - Response data automatically converted to pandas DataFrame format
        """
        assert self.url, "CODERS API URL is not configured. Please check your configuration."
        self.table_query= f"{self.url}/{table_name}{self.api_query}"
        self.table_response = requests.get(self.table_query)
        
        if self.table_response.status_code == 200:
            return pd.DataFrame.from_dict(self.table_response.json())
        else:
            raise RuntimeError(f">> Error fetching data for {table_name}: {self.table_response.status_code}")

    
    
    
    def load_local_data(self, 
                        table_name: str, 
                        region_code: str = None) -> pd.DataFrame:
        
        """
        Load data from a local file if it exists.
        
        ### Args:
            If region set to NONE, loads Country data.
        """

        file_name = f"{table_name}.pkl" if region_code is None else f"{table_name}_{region_code}.pkl"
        file_path = Path(self.data_pull['root']) / self.data_pull.get(table_name)/file_name
        file_path.parent.mkdir(parents=True, exist_ok=True)  # Creates parent directories if not exists.

        if file_path.is_file():
            print(f">> Loading data from local file: {file_path}")
            return pd.read_pickle(file_path)
        else:
            print(f">> No local file found at: {file_path}")
            return None  # Return None if the file does not exist

    def save_data(self, 
                  data: pd.DataFrame|gpd.GeoDataFrame, 
                  table_name: str, 
                  region_code: str = None):
        
        """Save the fetched data to a pkl file."""
        file_name = f"{table_name}.pkl" if region_code is None else f"{table_name}_{region_code}.pkl"
        file_path = Path(self.data_pull['root']) / self.data_pull.get(table_name)/file_name
        
        data.to_pickle(file_path)
        print(f"{table_name} data saved to:\n {file_path}")

    def create_gdf(self, df: pd.DataFrame) -> gpd.GeoDataFrame:
        """Create a GeoDataFrame from the given DataFrame."""
        df = df.copy()
        # Create a geometry column
        df['geometry'] = df.apply(lambda row: Point(row['longitude'], row['latitude']), axis=1)

        # Convert the DataFrame to a GeoDataFrame
        gdf = gpd.GeoDataFrame(df, geometry='geometry', crs='EPSG:4326')
        return gdf

    def get_table_canada(self, table_name: str, force_update: bool = False):
        """Get generator data for all of Canada.
        
        Args:
            table_name (str): The name of the table to fetch data from.
            force_update (bool): If True, force a data fetch from the API, ignoring local data.

        Returns:
            Tuple[pd.DataFrame, gpd.GeoDataFrame]: The generator data as a DataFrame and as a GeoDataFrame.
        """
                   
        file_path = Path(self.data_pull['root']) / self.data_pull.get(table_name, f"{table_name}.pkl")
        file_path.parent.mkdir(parents=True, exist_ok=True)  # Creates parent directories if not exists.
        
        # Check if the data file exists locally and if force_update is not set
        if file_path.is_file() and not force_update:
            data = pd.read_pickle(file_path)  # Load from local CSV
            print(f"Loaded {table_name} data from local file: {file_path}")
        else:
            # Fetch data from API if not found locally or if force_update is set
            data = self.fetch_data(table_name)
            print(f">> Data pulled {table_name} from [source checked: CODERS(https://sesit.dev/api/docs)]")
            self.save_data(data, table_name)
            
        df=data
        
        # Check if table_name contains "lines"; if it does, skip creating the GeoDataFrame
        if "lines" not in table_name:
            gdf = self.create_gdf(data)  # Only create GeoDataFrame if "lines" is not in table_name
        else:
            gdf = gpd.GeoDataFrame()  # Or however you wish to handle this case
        return df,gdf


    def get_table_provincial(self, 
                             table_name,
                             region_code: str = 'BC',
                             force_update: bool = False):
        """Get generator data for a specific province.
        
        Args:
            table_name (str): The name of the table to fetch data from e.g. 'substations','transmission_lines','generators' etc.
            force_update (bool): If True, force a data fetch from the API, ignoring local data.
        """
        df,gdf= self.get_table_canada(table_name, force_update=force_update)
            

        if "lines" not in table_name:
            # Apply provincial mask
            if "lines" not in table_name:
                data = gdf
            else:
                data = df

            region_mask = data['province'] == region_code
            self.region_data = data[region_mask]

            if not self.region_data.empty:
                return self.region_data  # Return the filtered GeoDataFrame
            else:
                return None


    def pull_all_data(self):
        
        for source, tables in self.data_pull.items():
            for table in tables:
                # Fetch the data using your existing methods
                data = self.fetch_data(table, source)
                print(f">> Retrieved data for {source}/{table} and saving to store")
                
                # Store the fetched data
                self.store.to_store(data, f'{source}/{table}',force_update=True)
         
    @staticmethod
    def convert_to_geodata(df):
        df['geometry'] = df.apply(lambda row: Point(row['longitude'], row['latitude']), axis=1)
        return gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326")

    def save_scenario_data(self, data: pd.DataFrame, table_name: str, scenario_col: str, folder: str):
        for scenario in data[scenario_col].unique():
            scenario_data = data[data[scenario_col] == scenario]
            file_path = self.data_pull_root / folder / f"{table_name}_{scenario}.csv"
            file_path.parent.mkdir(parents=True, exist_ok=True)
            scenario_data.to_csv(file_path, index=False)
            print(f"{table_name} scenario data saved to: {file_path}")

    def save_visualization(self, df: pd.DataFrame, x_col: str, y_cols: list, title: str, output_folder: str):
        fig = px.area(df, x=x_col, y=y_cols, title=title, labels={x_col: 'Year', 'value': 'Pj'})
        file_path = Path(output_folder) / f"{self.CRC}_visualization.html"
        file_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(file_path))
        print(f"Visualization saved to: {file_path}")