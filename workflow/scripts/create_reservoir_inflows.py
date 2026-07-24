from Z_legacy.pypsa_bc import hydro
from Z_legacy.pypsa_bc import utils
import atlite
import pandas as pd
import geopandas as gpd

# handles the config loading centrally
from Z_legacy.pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()

def main():
    # Description: main script for creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.

    # (i) get configuration  
    cfg = pypsa_aparser.pypsa_cfg

    utils.print_update(level=1,message="Preparing inflows for reservoirs...")


    # (ii) Read basin and site data. (Basins NA and Artic)
    na_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["na_file"])
    ar_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["artic_file"])
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))
    cutout = atlite.Cutout(path=utils.get_cutout_path(cfg))

    # 1) Load in hydroelectric generation sites
    # 2) Load in reservoirs
    hydro_sites = hydro.load_hydro_sites(cfg["output"]["create_hydro_assets"]["hydro_generation"])
    reservoir_sites = hydro.load_reservoir_sites(cfg["output"]["create_hydro_assets"]["hydro_reservoir"])

    # Create and save reservoir inflows
    hydro.create_cascade_inflow(reservoir_sites, basin_data, cutout,
                    hydro_sites, cfg, 
                    method = cfg["output"]["reservoir_inflows"]["inflow_method"])
    

    
if __name__ == '__main__':
    main()