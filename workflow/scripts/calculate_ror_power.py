from bc_power import hydro
from bc_power import utils
import atlite
import pandas as pd
import geopandas as gpd
from collections import namedtuple




def main():
    # Description: main script for creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.
    # NOTE: This is an old script and not in use anymore
    # (i) get configuration
    config_file = r"config/config2.yaml"
    cfg = utils.load_config(config_file)

    # (i) read basin and site data
    # Basins NA and artica
    na_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["na_file"])
    ar_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["artic_file"])
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))
    cutout = atlite.Cutout(path=cfg["cutout"]["file"])

    # Sites/plants which need inflow data. Reduce assets to only asset_id since asset_id is 1-to-1 with inflows
    sites_prep = hydro.load_hydro_sites(cfg["ror_inflows"]["ror_assets"])
    # hydro.rename_duplicate_asset_id(sites_prep)
    subset = ["asset_id","latitude","longitude",] 
    sum_list = ["capacity", "annual_avg_energy", "ramp_up", "ramp_down"] 
    sites = sites_prep.groupby(by="asset_id", group_keys=False).apply(lambda x: hydro.merge_assets(x, subset, sum_list))

    # (ii) Calculated the inflows for each site
    basins = hydro.prepare_basins(sites, basin_data)
    basin_inflows = hydro.calculate_basin_inflows(basins, cutout, height=bool(cfg['ror_inflows']['height']))
    site_inflows = hydro.calculate_plant_inflows(basin_inflows, basins, flowspeed=cfg['ror_inflows']['flowspeed'])

    # (iii) Compute the power availability series for each site based on capacity, inflow series, and
    # annual energy production
    power_series = hydro.calculate_ror_power(sites, site_inflows)
    power_series.to_csv(cfg['ror_inflows']['ror_outfile'])


if __name__ == '__main__':
    main()