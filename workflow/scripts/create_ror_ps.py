from pypsa_bc import hydro
from pypsa_bc import utils
import atlite
import pandas as pd
import geopandas as gpd
import sys
from collections import namedtuple
from pathlib import Path

def main(config_file:str|Path):
    # Description: main script for creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.

    # (i) get configuration
    # config_file = r"config/config.yaml"   
    cfg = utils.load_config(config_file)

    utils.print_update(level=1,message="Preparing existing ROR assets...")
  

    # (ii) Read basin and site data. (Basins NA and Artic)
    na_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["na_file"])
    ar_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["artic_file"])
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))
    cutout = atlite.Cutout(path=utils.get_cutout_path(cfg))

    # 1) Load in hydroelectric generation sites
    hydro_sites = hydro.load_hydro_sites(cfg['output']["create_hydro_assets"]["hydro_generation"])

    #NOTE: Alberta Special remove irrigation hydroelectric sites
    # irrigation_site_ids = ['AB_CHI_GSS','AB_ICP_GSS','AB_RYM_GSS']
    # hydro_sites = hydro_sites[hydro_sites['asset_id'].apply(lambda x: x not in irrigation_site_ids)].reset_index(drop=True)
    # 1) Create RoR power availability series:
    all_ror_sites = hydro_sites[(hydro_sites["hydro_type"] == "ror") | (hydro_sites["hydro_type"] == "ror-water")]
    hydro.create_ror_power(all_ror_sites, basin_data, cutout, cfg)
    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python create_hydro_assets.py <config_file>")
        sys.exit(1)
    config_file = sys.argv[1]
    main(config_file)