from bc_power import hydro, utils
import atlite
import pandas as pd
import geopandas as gpd
from collections import namedtuple


def main():
    # Description: main script for creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.

    # (i) get configuration
    config_file = r"config/config2.yaml"
    cfg = utils.load_config(config_file)


    # (ii) Read basin and site data. (Basins NA and Artic)
    na_basin_data = hydro.load_hydro_basins(cfg["data"]["basin_files"]["na_file"])
    ar_basin_data = hydro.load_hydro_basins(cfg["data"]["basin_files"]["artic_file"])
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