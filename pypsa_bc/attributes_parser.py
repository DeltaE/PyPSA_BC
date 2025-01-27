import yaml
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict
import logging as log

# Local Package
from pypsa_bc import utils

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
    pypsa_config_path:str|Path= field(default='config/data.yaml')
    
    def __post_init__(self):

        self.pypsa_config_path=Path(self.pypsa_config_path)

        # Define the path and filename
        # self.store = Path("data/store/BC_Combined_Modelling_results_.h5")
        # self.store.parent.mkdir(parents=True, exist_ok=True)
        
        self.data_cfg_path:str|Path=('config/data.yaml')
        self.pypsa_cfg:dict=self.load_config(self.pypsa_config_path)
        self.log = log.getLogger(__name__)
        self.check_dirs
    
    @staticmethod
    def load_config(config_file_path)->dict:
        """ 
        Loads the yaml file as dictionary and extracts the attributes to pass on child classes. 
        """
        with open(config_file_path, 'r') as file:
            data:dict = yaml.safe_load(file)
        return data
    
    @property
    def get_snapshot(self)->tuple:
        """
        Extracts the snapshot from the User Config file.
        """
        start_time = self.pypsa_cfg['province_mapping']['BC']['snapshots_tz_BC']['start'][0]
        end_time = self.pypsa_cfg['province_mapping']['BC']['snapshots_tz_BC']['end'][0]
        return (start_time, end_time)
    
    @property
    def get_visual_root(self):
        visual_root=Path('vis/pypsa')
        visual_root.mkdir(parents=True, exist_ok=True)
        return visual_root
    
    @property
    def check_dirs(self):
        directories = [
            "data/processed_data/load",
            "data/processed_data/network",
            "data/processed_data/wind",
            "data/processed_data/hydro",
            "data/processed_data/solar",
            "data/processed_data/tpp",
            "data/pypsa_data",
            'vis',
            'results/pypsa',
        ]
        
        for directory in directories:
            path = Path(directory)
            path.mkdir(parents=True, exist_ok=True)