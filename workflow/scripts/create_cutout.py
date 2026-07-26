from pathlib import Path
from zipfile import ZipFile

import geopandas as gpd
import pandas as pd
import pygadm
from pypsa_bc import attributes_parser, hydro, utils

aparser = attributes_parser.AttributesParser()


def ensure_shp_from_zip(zip_path: Path, shp_path: Path, remove_zip: bool = True) -> Path:
    """Extract one shapefile bundle from zip into the zip directory and return the .shp path."""
    zip_path = Path(zip_path)
    shp_path = Path(shp_path)

    # If shp already exists, prefer it and avoid touching zip.
    if shp_path.exists():
        return shp_path

    if not zip_path.exists():
        raise FileNotFoundError(f"Zip file not found: {zip_path}")

    out_dir = zip_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    with ZipFile(zip_path, "r") as zf:
        shp_members = [n for n in zf.namelist() if n.lower().endswith(".shp")]
        if not shp_members:
            raise FileNotFoundError(f"No .shp found inside zip: {zip_path}")

        shp_member = None
        for member in shp_members:
            if Path(member).name == shp_path.name:
                shp_member = member
                break
        if shp_member is None:
            shp_member = shp_members[0]

        base = Path(shp_member).stem
        parent = Path(shp_member).parent.as_posix()

        if parent in ("", "."):
            related = [
                n for n in zf.namelist()
                if Path(n).stem == base
            ]
        else:
            related = [
                n for n in zf.namelist()
                if Path(n).parent.as_posix() == parent and Path(n).stem == base
            ]

        # Write only same-stem files directly beside the zip.
        for member in related:
            target = out_dir / Path(member).name
            with zf.open(member) as src, target.open("wb") as dst:
                dst.write(src.read())

    extracted_shp = out_dir / Path(shp_member).name
    if not extracted_shp.exists():
        raise FileNotFoundError(f"Extracted .shp not found at expected path: {extracted_shp}")

    if remove_zip and zip_path.exists():
        zip_path.unlink()

    return extracted_shp

def get_basin_data(data_cfg):
    # from workflow.scripts import fetch_inputs
    from . import fetch_inputs  # Importing fetch_inputs directly to avoid circular import issues
    
    na_basin_zip_path = Path(data_cfg["remote"]["HydroBASINS"]['dest'])
    ar_basin_zip_path = Path(data_cfg["remote"]["HydroBASINS_arctic"]['dest'])
    na_basin_shp_path = Path(data_cfg["basin_files"]["na_file"])
    ar_basin_shp_path = Path(data_cfg["basin_files"]["arctic_file"])

    if not na_basin_shp_path.exists() and not na_basin_zip_path.exists():
        fetch_inputs.main(only=["HydroBASINS"])
    if not ar_basin_shp_path.exists() and not ar_basin_zip_path.exists():
        fetch_inputs.main(only=["HydroBASINS_arctic"])
    
    ensure_shp_from_zip(na_basin_zip_path, na_basin_shp_path, remove_zip=True)
    ensure_shp_from_zip(ar_basin_zip_path, ar_basin_shp_path, remove_zip=True)

    na_basin_data = hydro.load_hydro_basins(na_basin_shp_path) 
    ar_basin_data = hydro.load_hydro_basins(ar_basin_shp_path)
    
    basin_data = gpd.GeoDataFrame(pd.concat([na_basin_data, ar_basin_data]))

    return basin_data

def get_country_boundaries(data_cfg:dict, level:int=1):
    gadm_file = Path(data_cfg["GADM"][f"country_file_L{level}"])
    gadm_file.parent.mkdir(parents=True, exist_ok=True)

    if gadm_file.exists():
        return gpd.read_file(gadm_file)

    gdf = pygadm.Items(name="Canada", content_level=level)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")  # CRS degrees
    gdf.to_file(gadm_file, driver="GeoJSON")
    return gdf

def audit_cutout(cfg) -> pd.DataFrame:
        """
        Inspect if a cutout file already exists on disk for the specified year.
        Restores the instance's original snapshot state after probing.
        """
        _org_path=str(utils.get_cutout_path(cfg))
        _org_year=str(cfg['cutout']["snapshots"]["start"][0][:4])
        path_template = _org_path.replace(_org_year, '{year}')

        rows = []
        for year in range(int(cfg['cutout']["snapshots"]["start"][0][:4]),
                          int(cfg['cutout']["snapshots"]["end"][0][:4]) + 1):
            p       = Path(path_template.format(year=year))
            exists  = p.exists()
            size_mb = round(p.stat().st_size / 1e6, 1) if exists else 0
            rows.append({"year": year, "path": str(p), "exists": exists, "size_MB": size_mb})


        return pd.DataFrame(rows)
    
def main(): 
    
    # Description: main script for creating a cutout of the modelled region
    # creating the hydro cutout based on hydro site locations and the basins they are located within
    # and each basins upstream basins.
    cfg =aparser.data_cfg
    cutout_cfg=aparser.params_cfg
    
    print("Preparing cutout...")
    
    basin_data = get_basin_data(cfg)
    
    hydro_sites = hydro.load_hydro_sites(cfg["output"]["create_hydro_assets"]["hydro_generation"])
    # Get bounds needed for hydro cutout: Polygon containing smallest bounds to contain all hydro basins which 
    # provide water to the hydro assets of interest.
    hydro_polygon = hydro.get_hydro_cutout_polygon(hydro_sites, basin_data)

    # Determine true largest bounds based on max/min of hydro_bounds and the regional bounds
    gdf = get_country_boundaries(cfg,level=1)
    
    # region_map_2_gadm = {"BC":"BritishColumbia", "AB":"Alberta",
    #                      "SK":"Saskatchewan", "MB":"Manitoba"}  # update needed with region mapping set in cfg
    
    region_list = ['CA-BC']#[region_map_2_gadm[region] for region in cfg['output']['prepare_base_network']['regions']]
    
    mask = gdf['ISO_1'].apply(lambda x: x in region_list)
    region_polygon = utils.get_region_polygon(gdf[mask].geometry)

    # get max/min bounds to create a bounding box for the cutout
    bounds = utils.get_bounds([hydro_polygon, region_polygon])

    # Get resolution for ERA5
    if cutout_cfg["cutout"]["module"][0] == "era5":
        utils.create_era5_cutout(bounds, cutout_cfg)
    else:
        source = cutout_cfg["cutout"]["module"]
        print(f"Creating cutouts for {source} has not been implemented yet!")

if __name__ == '__main__':
    main()
    