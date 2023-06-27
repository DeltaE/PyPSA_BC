from bc_power import hydro
from bc_power import utils
import pandas as pd
import geopandas as gpd
import atlite

def main():
    # Description: main script for creating a cutout of the modelled region
    # creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.
    # Configuration inputs:
    config_file = r"/mnt/c/Users/pmcw9/Delta-E/PICS/PyPSA_BC/config/config.yaml"

    # main():
    # (i) read needed data
    cfg = utils.load_config(config_file)

    # By default read in hydro data
    na_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["na_file"])
    ar_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["artic_file"])
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))
    hydro_sites = hydro.load_hydro_sites(cfg["hydro_prep"]["hydro_generation"])
    # Get bounds needs for hydro cutout
    hydro_polygon = hydro.get_hydro_cutout_polygon(hydro_sites, basin_data)

    # Determine true largest bounds based on max/min of hydro_bounds and the regional bounds
    gdf = gpd.read_file("/mnt/c/Users/pmcw9/Delta-E/PICS/Data/regions/gadm41_CAN_1.json")
    mask = gdf['NAME_1'] == "BritishColumbia"
    region_polygon = utils.get_region_polygon(gdf[mask].geometry)

    # get max/min bounds to create a bounding box for the cutout
    bounds = utils.get_bounds([hydro_polygon, region_polygon])

    # Get resolution for ERA5
    if cfg["cutout"]["source"] == "era5":
        utils.create_era5_cutout(bounds, cfg)
    else:
        source = cfg["cutout"]["source"]
        print(f"Creating cutouts for {source} has not been implemented yet!")

if __name__ == '__main__':
    main()