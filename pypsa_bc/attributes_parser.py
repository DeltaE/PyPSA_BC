import yaml
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict
import logging as log

# Local Package
from .clews import model_structure as clews_const
from . import utils

# Logging Configuration
log.basicConfig(level=log.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

"""
# Key Changes and Benefits over v1
    The @dataclass decorator simplifies class creation and automatically generates the __init__, __repr__, and other methods.

## Field Initialization:
    Attributes that require processing during initialization (like reading configurations) are defined with init=False and processed in the __post_init__ method.

## Default Values:
    The resource_type has a default value specified directly in the field declaration, which simplifies the __init__ method.

## Type Annotations:
    Type hints enhance code readability and help with type checking tools.
"""
# Logging Configuration
log.basicConfig(level=log.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass
class AttributesParser:
    """
    This is the parent class that will extract the core attributes from the User Config file.
    
    ## Remarks:
        - Handles the clews_builder.yaml file creation and updates under-the hood.  
        - The clews_builder.yaml file if not exists will mirror the skeleton file sourced from 'models/BC_Nexus/config/clews_builder_skeleton.yaml'
    """
    # Attributes that are required as Args.
    combined_model_config_path:str|Path= field(default='config/config.yaml')
    
    def __post_init__(self):

        self.combined_model_config_path=Path(self.combined_model_config_path)

        # Define the path and filename
        # self.store = Path("data/store/BC_Combined_Modelling_results_.h5")
        # self.store.parent.mkdir(parents=True, exist_ok=True)
        
        self.data_cfg_path:str|Path=('config/data.yaml')
        self.clews_builder_config_path:str|Path=Path('config/clews_builder.yaml') # Kept it as under the hood config file build and updated automatically.
        self.clews_builder_skeleton_source:str|Path=Path('models/BC_Nexus/config/clews_builder_skeleton.yaml')
        self.build_clews_builder_skeleton()
        
        # Load the user configuration master file by using the method
        self.data_cfg:Dict[str:dict]=AttributesParser.load_config(self.data_cfg_path)
        self.cm_config:Dict[str,dict] = AttributesParser.load_config(self.combined_model_config_path)
        self.clewsb_config:Dict[str,dict] = AttributesParser.load_config(self.clews_builder_config_path)
        self.clews_cm:Dict= self.cm_config.get('clews','clews')
        
        self.log = log.getLogger(__name__)
    
    def build_clews_builder_skeleton(self):
        if not self.clews_builder_config_path.exists():
            self.clews_builder_skeleton: dict = AttributesParser.load_config(self.clews_builder_skeleton_source)
            with self.clews_builder_config_path.open('w') as file:
                yaml.dump(self.clews_builder_skeleton, file)
                utils.print_update(level=2,message=f"clews builder config was missing, created at {self.clews_builder_config_path} with the skeleton.")
        else:
            pass
    
    def get_LandCluster_paths(self)->tuple:
        source='models/BC_Nexus/LandClusterData'
        dest='data/clews_data/LandClusterData'
        return (source,dest)
    
    def get_SETs_Path(self):
        return Path('data/clews_data/SETs')
    
    @staticmethod
    def load_config(config_file_path):
        """ 
        Loads the yaml file as dictionary and extracts the attributes to pass on child classes. 
        """
        with open(config_file_path, 'r') as file:
            data = yaml.safe_load(file)
        return data
    
    def get_SETs_save_to(self):
        SETS_save_to=Path('data/clews/SETs')
        SETS_save_to.mkdir(parents=True,exist_ok=True)
        return SETS_save_to
    
    
    def get_clews_snapshot(self) -> tuple[int, int]:
        # start_year = int(self.clews_cm['GENERAL'].get('start_year', 2021))  # Default to 2021 if missing
        # last_year = int(self.clews_cm['GENERAL'].get('last_year', 2050))  # Default to 2050 if missing
        clews_snapshot=clews_const.snapshot
        start_year=clews_snapshot['start']
        last_year=clews_snapshot['end']
        return (start_year, last_year)

    def get_region(self):
        # region=self.clews_cm['GENERAL'].get('region','REGION1')
        clews_regions:dict=clews_const.Regions
        current_region=next(iter(clews_regions))
        return current_region
    
    def get_default_csv_template_path(self)->Path:
        # csv_template_path=Path(self.clewsb_config.get('clews_datapackage_template','models/BC_Nexus/csv_template'))
        csv_template_path=Path('data/clews_data/csv_template')
        
        return csv_template_path
    
    def get_input_data_csv_dir(self)->Path:
        input_csv_dir=Path('data/clews_data/inputs_csv')
        input_csv_dir.mkdir(parents=True,exist_ok=True)
        return input_csv_dir
    
    def get_otoole_yaml_file(self,
                             storage_algorithm:str='Kotzur')->Path:
        
        return Path(f'models/BC_Nexus/Model_{storage_algorithm}/otoole_config_{storage_algorithm}.yaml')
    
    def get_linking_tool_results(self):
        return Path('results/linking')
        
    def get_scenarios(self,
                      clews_cm:dict=None)->dict:
        return clews_cm['SCENARIOS'] if clews_cm else self.clews_cm['SCENARIOS']
        
    def get_data_path_attributes(self):
        """
        ## Returns:
        
            - tuple: A tuple containing the following elements in order:  

            0 → `StorageMaxCapacity` (dict)  
                - Maximum storage capacity available for different regions or technologies.  

            1 → `ResidualStorageCapacity` (dict)  
                - Remaining or existing storage capacity after accounting for prior allocations.  

            2 → `ext_wind_CF` (pd.DataFrame)  
                - External wind capacity factors, providing wind resource availability over time.  

            3 → `ext_solar_CF` (pd.DataFrame)  
                - External solar capacity factors, representing solar generation potential over time.  

            4 → `demand_profile` (pd.DataFrame)  
                - Time-series data of energy demand across different regions or sectors.  

            5 → `hydro_res_CF` (pd.DataFrame)  
                - Capacity factors for hydro reservoir plants, indicating generation potential.  

            6 → `hydro_ror_ts` (pd.DataFrame)  
                - Time-series data for run-of-river hydro generation availability.  

            7 → `future_wind_CF` (pd.DataFrame)  
                - Projected wind capacity factors for future resource options.  

            8 → `future_solar_CF` (pd.DataFrame)  
                - Projected solar capacity factors for future resource options.  

            9 → `committed_wind_CF` (pd.DataFrame)  
                - Wind capacity factors for already committed wind projects.  

            10 → `committed_solar_CF` (pd.DataFrame)  
                - Solar capacity factors for already committed solar projects.  

            11 → `clews_CF_csv_file` (str)  
                - Path to the CLEWs capacity factor CSV file containing processed capacity factor data.  

            12 → `clews_demand_profile` (str)  
                - Path to the CLEWs demand profile CSV file with processed demand data. ` 
        """

        linking_tool_results=self.get_linking_tool_results()
        
        StorageMaxCapacity =  self.clewsb_config['STORAGE']['BATTERY']['storage_max_capacity']
        ResidualStorageCapacity =  self.clewsb_config['STORAGE']['BATTERY']['residual_storage_capacity']
        
        # Existing Resources' Data
        ext_wind_CF = self.data_cfg['pypsa']['output']['create_ext_wind_ts']['fname'] # Path('data/processed_data/wind/existing/bc_ext_wind_ts.csv') 
        ext_solar_CF = self.data_cfg['pypsa']['output']['create_ext_solar_ts']['fname'] # Path('data/processed_data/solar/existing/bc_ext_solar_ts.csv') 
        demand_profile =  Path('data/downloaded_data/CODERS/data-pull/demand/BC_provincial_demand_profile.csv') 
        hydro_res_CF = self.data_cfg['pypsa']['output']['reservoir_inflows']['fname']
        hydro_ror_ts=self.data_cfg['pypsa']['output']['ror_ps']['fname']
        
        # Future Resource Options (Committed/Planned)
        future_wind_CF=linking_tool_results / 'resource_options_wind_timeseries.csv'
        future_solar_CF=linking_tool_results /'resource_options_solar_timeseries.csv'
        
        # Committed Resource Options, contracted but not existing (i.e. we don't have existing profiles for them)
        committed_wind_CF=Path('results/linking/BCH_CFP24_wind_ts.csv')
        committed_solar_CF=Path('results/linking/BCH_CFP24_solar_ts.csv')
        
        # Profiles (Capacity Factort, Specific Demand) from Existing Resources, formatted in CLEWs schema
        clews_CF_csv_file =  Path('data/clews_data/inputs_csv_8760/CapacityFactor.csv')
        clews_CF_csv_file.parent.mkdir(parents=True, exist_ok=True)
        
        clews_demand_profile =  Path('data/clews_data/inputs_csv_8760/SpecifiedDemandProfile.csv') 
        clews_demand_profile.parent.mkdir(parents=True, exist_ok=True)
        
        return (
            StorageMaxCapacity,
            ResidualStorageCapacity,
            ext_wind_CF,
            ext_solar_CF,
            demand_profile,
            hydro_res_CF,
            hydro_ror_ts,
            future_wind_CF,
            future_solar_CF,
            committed_wind_CF,
            committed_solar_CF,
            clews_CF_csv_file,
            clews_demand_profile       
        )
    
    def get_case_input_csvs_path(
        self, 
        storage_algorithm: str
    ):
        
        case_input_csvs_path= Path(f'data/clews_data/Model_{storage_algorithm}/inputs_csv')
        case_input_csvs_path.mkdir(parents=True,exist_ok=True)
        
        return case_input_csvs_path
    
    def get_scenario_results_path(self, 
                                  scenario:str,
                                  storage_algorithm:str,):
        
        scenario_results_root:Path=Path(f'results/clews/Model_{storage_algorithm}_{scenario}')
        scenario_results_root.mkdir(parents=True,exist_ok=True)
        
        return scenario_results_root
    
    def get_scenario_inputs_path(self,
                                scenario:str,
                                storage_algorithm:str,
                               ):
        
        case_ip_csvs=self.get_case_input_csvs_path(storage_algorithm)
        
        scenario_inputs_root:Path=case_ip_csvs.parent/f'{scenario}'
        scenario_inputs_root.mkdir(parents=True,exist_ok=True)

        return  scenario_inputs_root
        
    def get_model_file_path(self,
                            storage_algorithm:str):
        
        org_model_file = Path(f'models/BC_Nexus/Model_{storage_algorithm}/model_BCNexus_{storage_algorithm}.txt')

        
        return org_model_file
    """,
               storage_case_model_file_preprocessed,
               ip_data_file,
               storage_case_ip_data_file_preprocessed)
               """
           
    def get_data_and_model_file_paths(self,
                                      input_csvs:str|Path,
                                      storage_algorithm:str,
                                      scenario:str):
        
        org_model_file = Path(f'models/BC_Nexus/Model_{storage_algorithm}/model_BCNexus_{storage_algorithm}.txt')
        # ip_data_file= input_csvs.parent/f'BCNexus_{storage_algorithm}.txt' 
        
        # (scenario_inputs_root,
        # scenario_results_root,
        # LP_file) = self.get_scenario_attributes(scenario,storage_algorithm)
        
        # storage_case_model_file_preprocessed = scenario_inputs_root/f'model_BCNexus_{storage_algorithm}_preprocessed.txt'
        # storage_case_ip_data_file_preprocessed = scenario_inputs_root/f'BCNexus_{storage_algorithm}_processed.txt'
        
        return(org_model_file)
    """,
               storage_case_model_file_preprocessed,
               ip_data_file,
               storage_case_ip_data_file_preprocessed)
               """