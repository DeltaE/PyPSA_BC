from bc_power import hydro
from bc_power import utils
import atlite
import pandas as pd
import geopandas as gpd
from collections import namedtuple




def main():
    # Description: main script for creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.

    # (i) get configuration
    config_file = r"/mnt/c/Users/pmcw9/Delta-E/PICS/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)


    # (ii) Read basin and site data. (Basins NA and Artic)
    na_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["na_file"])
    ar_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["artic_file"])
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))
    cutout = atlite.Cutout(path=cfg["hydro_prep"]["cutout"])

    # 1) Load in hydroelectric generation sites
    # 2) Load in reservoirs
    hydro_sites = hydro.load_hydro_sites(cfg["reservoir_inflows"]["hydro_assets"])
    reservoir_sites = hydro.load_reservoir_sites(cfg["reservoir_inflows"]["reservoir_assets"])

    # (iii) Determine which reservoirs are 
    # a) Are apart of the cascaded reservoirs
    # b) Have WUP data
    # c) Head reservoir (i.e. no reservoirs up stream) OR downstream

    # a) Split into cascade and non-cascade
    rid_list = reservoir_sites.index.tolist()
    mask_cascade = hydro_sites['upper_reservoir_id'].apply(lambda x: x in rid_list)
    cascade_sites = hydro_sites[mask_cascade].copy()
    ror_sites = hydro_sites[~mask_cascade].copy()

    # b) Get mask of cascades reservoirs with WUP statistics
    mask_res_wup = cascade_sites['upper_reservoir_id'].apply(lambda x: hydro.check_wup_exists(x,
                                     cfg['bc_hydro']['inflow_tables']))

    # c) Items in cascade_sites[~mask_res_wup] need to be checked for upstream or downstream
    # upstream head reservoirs (no reservoirs upstream of these) with no WUP statistics
    mask_res_up_no_wup = cascade_sites['upper_reservoir_id'].apply(lambda x: 
                        hydro.check_head_reservoir(x, cascade_sites)) & (~mask_res_wup)
    # downstream reservoirs with no WUP statistics
    mask_res_down_no_wup = (~cascade_sites['upper_reservoir_id'].apply(lambda x: 
                        hydro.check_head_reservoir(x, cascade_sites))) & (~mask_res_wup)

    # (iv) Save the inflow series for reservoirs and power availability series for generation assets
    # 1) Reservoir inflow normalized with WUP stats: cascade_site[mask_res_wup]
    # 2) Rerservoir inflow set to 0: Some of cascade_sites[mask_res_down_no_wup]
    rid_up_down_no_wup = cascade_sites[mask_res_down_no_wup]['upper_reservoir_id'].unique().tolist()
    rid_res_wup =  cascade_sites[mask_res_wup]['upper_reservoir_id'].unique().tolist()
    rid_2_method = {'Normalize':rid_res_wup, "Impute":rid_up_down_no_wup}

    hydro.create_cascade_inflow(reservoir_sites, basin_data, cutout,
        cascade_sites, rid_2_method, cfg, method = "mean_inflow_normalize")

    # 3) Create RoR power availability series: cascade_sites[mask_res_up_no_wup] AND ror_sites
    all_ror_sites = gpd.GeoDataFrame(pd.concat([ror_sites, cascade_sites[mask_res_up_no_wup]]))
    hydro.create_ror_power(all_ror_sites, basin_data, cutout, cfg)


if __name__ == '__main__':
    main()