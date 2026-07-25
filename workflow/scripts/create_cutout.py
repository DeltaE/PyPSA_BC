from pypsa_bc import hydro
from pypsa_bc import utils
import pandas as pd
import geopandas as gpd


def main(): 
    
    # Description: main script for creating a cutout of the modelled region
    
    # creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.
    # Configuration inputs:
    config_file = r"config/data.yaml"
    cfg_cmplt = utils.load_config(config_file)
    cfg=cfg_cmplt['pypsa']
    
    print("Preparing cutout...")
    
    # By default read in hydro data
    na_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["na_file"]) 
    ar_basin_data = hydro.load_hydro_basins(cfg["basin_files"]["artic_file"])
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))
    hydro_sites = hydro.load_hydro_sites(cfg["output"]["create_hydro_assets"]["hydro_generation"])
    # Get bounds needed for hydro cutout: Polygon containing smallest bounds to contain all hydro basins which 
    # provide water to the hydro assets of interest.
    hydro_polygon = hydro.get_hydro_cutout_polygon(hydro_sites, basin_data)

    # Determine true largest bounds based on max/min of hydro_bounds and the regional bounds
    gdf = gpd.read_file(cfg_cmplt['GADM']['country_file_L1'])
    
    region_map_2_gadm = {"BC":"BritishColumbia", "AB":"Alberta",
                         "SK":"Saskatchewan", "MB":"Manitoba"}  # update needed with region mapping set in cfg
    
    region_list = [region_map_2_gadm[region] for region in cfg['output']['prepare_base_network']['regions']]
    
    mask = gdf['NAME_1'].apply(lambda x: x in region_list)
    region_polygon = utils.get_region_polygon(gdf[mask].geometry)

    # get max/min bounds to create a bounding box for the cutout
    bounds = utils.get_bounds([hydro_polygon, region_polygon])

    # Get resolution for ERA5
    if cfg_cmplt["cutout"]["module"][0] == "era5":
        utils.create_era5_cutout(bounds, cfg_cmplt)
    else:
        source = cfg_cmplt["cutout"]["module"]
        print(f"Creating cutouts for {source} has not been implemented yet!")

if __name__ == '__main__':
    main()
    