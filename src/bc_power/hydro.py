import geopandas as gpd
import atlite
from shapely.geometry import Point
import shapely
import pandas as pd
import numpy as np
from collections import namedtuple
import dask
from scipy.sparse import csr_matrix
import xarray as xr
from bc_power import utils
import os
from shapely.validation import make_valid


Basins = namedtuple("Basins", ["plants", "meta", "shapes"])

def load_hydro_basins(source_file):
    '''
    This function loads in a HydroBASINS shape file.
    source_file: Path + file name of the HydroBASINS file to be loaded.
    return: Returns a GeoDataFrame of the HydroBASINS shape file.
    '''
    return gpd.read_file(source_file).set_index('HYBAS_ID')

def load_hydro_sites(hydro_sites_file):
    '''
    This function loads in a prepared file of hydro sites.
    source_file: Path + file name of the hydro sites file to be loaded.
    return: Returns a GeoDataFrame of the hydro sites.
    '''
    gdf = gpd.read_file(hydro_sites_file).set_index('component_id')
    # Need to fix all at some point...
    gdf['latitude'] = gdf['latitude'].astype(float)
    gdf['longitude'] = gdf['longitude'].astype(float)
    gdf['capacity'] = gdf['capacity'].astype(float)
    gdf['annual_avg_energy'] = gdf['annual_avg_energy'].astype(float)
    gdf['cascade_order'] = pd.to_numeric(gdf['cascade_order'], errors='coerce')
    gdf.rename({'latitude':'lat', 'longitude':'lon'}, axis=1,inplace=True)

    return gdf

def load_reservoir_sites(hydro_sites_file):
    '''
    This function loads in a prepared file of hydro sites.
    source_file: Path + file name of the hydro sites file to be loaded.
    return: Returns a GeoDataFrame of the hydro sites.
    '''
    gdf = gpd.read_file(hydro_sites_file).set_index('asset_id')
    # Need to fix all at some point...
    gdf['latitude'] = pd.to_numeric(gdf['latitude'], errors='coerce')
    gdf['longitude'] = pd.to_numeric(gdf['longitude'], errors='coerce')
    gdf['min_storage'] = pd.to_numeric(gdf['min_storage'], errors='coerce')
    gdf['max_storage'] = pd.to_numeric(gdf['max_storage'], errors='coerce')
    gdf['min_level'] = pd.to_numeric(gdf['min_level'], errors='coerce')
    gdf['max_level'] = pd.to_numeric(gdf['max_level'], errors='coerce')
    gdf['cascade_order'] = pd.to_numeric(gdf['cascade_order'], errors='coerce')
    gdf.rename({'latitude':'lat', 'longitude':'lon'}, axis=1,inplace=True)

    return gdf

def get_site_basins(sites,meta,shapes):
    '''
    This functions determines which basin contains which hydroelectric site and all upstream basins.
    sites: GeoDataFrame of all the hydroelectric sites.
    meta: GeoDataFrame of all basin data except the geometry.
    shapes: GeoSeries taken from the geometry column of the basin data.
    '''
    
    plant_basins = [] 
    for site in sites.itertuples():
        # (i) find the basin using the hydrobasin shapes and the (lon,lat) of the site/plant
        # hid = hydro_index of the basin which intersects/contains the site/plant 
        # Assumes located in the first basin
        hid = shapes.index[shapes.intersects(Point(site.lon, site.lat))][0]

        # (ii) find the basins upstream (can always have multiple upstream)
        i = 0
        hids = [hid]
        while i < len(hids): # Finished when no more upstream...
            hids.extend(meta.index[meta["NEXT_DOWN"] == hids[i]])
            i += 1
            
        # (iii) Store all upstream hydrobasins
        plant_basins.append((hid,hids)) # (hydro_basin_index, all upstream hydro_basin_index)

    return pd.DataFrame(plant_basins, columns=["hid", "upstream"], index=sites.index)


def get_hydro_cutout_polygon(sites,basin_data):
    '''
    This function finds the min/max bounds needed to create a cutout for the hydrobasins.
    return: A shapely polygon, representing the smallest bounding box for the
            hydro basins this is then used later to find min/max bounds. 
    '''
    meta = basin_data[basin_data.columns.difference(("geometry",))]
    shapes = basin_data["geometry"]

    # Get DataFrame of each plant's basin and its corresponding upstream basins
    plant_basins = get_site_basins(sites, meta, shapes)
    unique_basins = pd.Index(plant_basins["upstream"].sum()).unique().rename("hid")
    basins = Basins(plant_basins, meta.loc[unique_basins], shapes.loc[unique_basins])

    # cutout bounds
    west_lon, south_lat, east_lon, north_lat = shapes.loc[basins.plants.hid[0]].bounds # initialize
    for ppl in basins.plants.itertuples():
        for hid in ppl.upstream:
            min_lon, min_lat, max_lon, max_lat = shapes.loc[hid].bounds
            west_lon = min(min_lon, west_lon)
            east_lon = max(max_lon, east_lon)
            south_lat = min(min_lat, south_lat)
            north_lat = max(max_lat, north_lat)

    bbox = (west_lon, south_lat, east_lon, north_lat)
    polygon = shapely.geometry.box(*bbox, ccw=True)
    
    return polygon

def create_cutout(cutout_file,bounds,dx,dy,time_horizon,features=['height','runoff'],module="era5"):
    '''
    This function creates a cutout from Atlite and saves it to the cutout_file (path + filename).
    NOT USED ANYMORE SINCE A SINGLE CUTOUT IS USED FOR EVERYTHING RATHER THAN CUTOUTS FOR EACH TECH.
    '''
    cutout = atlite.Cutout(path=cutout_file,
                        module=module,
                        x=slice(bounds['west_lon'] - dx, bounds['east_lon'] + dx),
                        y=slice(bounds['south_lat'] - dy,bounds['north_lat'] + dy ),
                        dx=dx,
                        dy=dy,
                        time=time_horizon,
                        prepared_features=features)

    cutout.prepare()

def merge_assets(df,subset,sum_list):
    '''
    Function used to reduce hydroelectric datasets from turbines to an aggregate asset.
    This aggregation operator is currently applied only to the installed capacities for the units and
    the annual_avg_energy.
    df: Dataframe which is passed via a groupby operation.
    subset: Name of columns to use for deduplication
    sum_list: Name of parameters/columns to aggegtate using the sum operation.
    Example:
    Input dataframe has following entries below:
    component_id | asset_id | capacity | annual_avg_energy
    BC_MCA01_GEN | BC_MCA_GSS | 492 | 1936.79
    BC_MCA02_GEN | BC_MCA_GSS | 492 | 1936.79
    BC_MCA03_GEN | BC_MCA_GSS | 494 | 1942.7
    BC_MCA04_GEN | BC_MCA_GSS | 494 | 1942.7
    BC_MCA05_GEN | BC_MCA_GSS | 500 | 1968.45
    BC_MCA06_GEN | BC_MCA_GSS | 500 | 1968.45

    Output dataframe will have the following:
    asset_id | capacity | annual_avg_energy
    BC_MCA_GSS | 2972 | 11695.88

    '''
    # Other columns don't matter here for calcualting inflow and associated power production.
    df_out = df.drop_duplicates(subset=subset).set_index("asset_id").copy()
    for param in sum_list:
        df_out[param] = df[param].sum()
    return df_out

def prepare_basins(sites, basin_data):
    '''
    This function will return an object linking the hid (basin ID from hydrBASINS) for the basin containing each sites/plant
    additionally it will include the name of each generator (component_id) to each 

    '''
    # seperate components of the basin_data
    meta = basin_data[basin_data.columns.difference(("geometry",))]
    shapes = basin_data["geometry"].apply(lambda x: make_valid(x))

    # get dataframe of the hid and upstream hids for each plant/site
    plant_basins = get_site_basins(sites, meta, shapes)

    # Unique basins indices of all upstream (inclusive of hid basin) basins
    # Cascade modification would happen here...
    unique_basins = pd.Index(plant_basins["upstream"].sum()).unique().rename("hid")

    # basins.plants, basins.meta, basins.shapes
    basins = Basins(plant_basins, meta.loc[unique_basins], shapes.loc[unique_basins])
    return basins


def prepare_cascade_basins(reservoir_sites, basin_data, gen_sites):
    '''
    Adjust upstream list of basins for each downstream reservoir within a cascade.
    This is done because when the upstream basins are found for reservoir downstream of cascades
    they will also include basins which are apart of the upstream reservoirs. This is not desired since
    inflow into these basins are going to be regulated by the upstream dams.
    reservoirs_sites: DataFrame of reservoirs sitings.
    basin_data: GeoDataFrame of all upstream basins.
    hydro_gen: Path + filename of the file containing the hydroelectric facilities.
               (This file indicates the ordering of cascades)
    '''
    # prepre base basins and load generation sites
    basins = prepare_basins(reservoir_sites, basin_data.copy())
    # gen_sites = load_hydro_sites(hydro_gen)

    # (i) create ordered dictionary based on stream relationship between site locations
    # i.e. {reservoir:set(direct upper reservoirs)}
    cascade_dict = {}
    for idx,row in gen_sites.iterrows():
        # error check for the empty lower 
        if row['lower_reservoir_id'] == "":
            continue
        if row['lower_reservoir_id'] not in cascade_dict.keys():
            cascade_dict[row['lower_reservoir_id']] = set([row['upper_reservoir_id']])
        else:
            cascade_dict[row['lower_reservoir_id']].add(row['upper_reservoir_id'])

    # (ii) Now add all the reservoir which are at the top
    for aid in reservoir_sites.index:
        if aid not in cascade_dict.keys():
            cascade_dict[aid] = set() # empty set since none are upstream of it

    # (iii) Now loop over reservoirs and modify the basins for each
    update_dict = {} # {aid:adjusted_hids}
    for aid,row in basins.plants.iterrows():
        if len(cascade_dict[aid]) == 0:
            continue
        else:
            # Create set of all upstream reservoir hids
            cascade_hids = set([])
            for up_aid in cascade_dict[aid]:
                cascade_hids.update(basins.plants.loc[up_aid,'upstream']) # inclusive

            # Determins hids between current reservoir and all other upstream reservoirs
            modified_hids = list(set(basins.plants.loc[aid,'upstream']) - cascade_hids)
            update_dict[aid] = modified_hids

    for aid,hids in update_dict.items():
        basins.plants.at[aid,'upstream'] = hids # new column otherwise

    return basins

def calculate_basin_inflows(basins, cutout, height=True):
    '''
    This function calculates the inflows for each basin. This is done by finding an average runoff for the each basin
    from the weighted average runoff from all interesecting cutouts. Averaging the proportion of runoff intersected by
    each basin by the cumulative proportions of cutouts intersected with the basin.
    More details can be found in a supplementary document.
    Units are returned in m^3/hr
    '''
    # index set before calling cutout.runoff()
    index = basins.shapes.index

    # Create indicators matrix:
    # columns represent the cutouts
    # rows represent the upstream basins
    # cell values indicates the ratio/proportion of cutout (j) area contained in the basin (i)
    matrix = cutout.indicatormatrix(basins.shapes)

    # normalize the matrix ...
    # normalizes the matrix based on the 
    # cell values indicates the proportion of intersected area between cutout (j) and basin (i)
    # w.r.t the total intersected areas of basin (i) with all cutouts (j)
    matrix_normalized = matrix / matrix.sum(axis=1)

    ## calculate runoff
    if height:
        # This is done for a potential energy approach for calculate power availability
        # In this case the inflow is not a true inflow.
        runoff = cutout.data['runoff'] * cutout.data['height'] 
    else:
        runoff = cutout.data['runoff'] 

    # Normalzied matrix...
    matrix_normalized = csr_matrix(matrix_normalized)

    # Aggregate the matrix
    if isinstance(runoff.data, dask.array.core.Array):
        da = runoff.stack(spatial=("y", "x")) # New dim which uses existing dims of (y,x)
        result = xr.apply_ufunc(
            lambda da: da * matrix_normalized.T,
            da,
            input_core_dims=[["spatial"]],
            output_core_dims=[[index.name]],
            dask="parallelized",
            output_dtypes=[da.dtype],
            dask_gufunc_kwargs=dict(output_sizes={index.name: index.size}),
        ).assign_coords(**{index.name: index})
    else:
        # need to update this error later
        print(f'ERROR: Runoff data is not in correct format, runoff data has {type(runoff.data)} type')
        exit(1)
        # da = runoff.stack(spatial=("y", "x")).transpose("spatial", "time")
        # result = xr.DataArray(matrix_normalized * da, [index, da.coords["time"]])
    result *= xr.DataArray(basins.shapes.to_crs(dict(proj="cea")).area)
    result.load()
    return result

def calculate_plant_inflows(basin_inflows, basins, flowspeed=1):
    '''
    Calculate inflow timeseries at hydroelectric plants by temporally shifting and aggregating
    the inflows at all upstream inflows for each basin.
    '''
    inflow = xr.DataArray(
            np.zeros((len(basins.plants), basin_inflows.indexes["time"].size)),
            [("plant", basins.plants.index), basin_inflows.coords["time"]],
        )

    for ppl in basins.plants.itertuples():
        # instantiate array for inflow to the plant
        inflow_plant = inflow.loc[dict(plant=ppl.Index)]
        # calculate distance from all upsteam basins to the current basin
        # DIST_MAIN: Units in km
        distances = (
            basins.meta.loc[ppl.upstream, "DIST_MAIN"]
            - basins.meta.at[ppl.hid, "DIST_MAIN"]
        )
        # Estimates time delay for inflow to reach the basin containing the power plant 
        nhours = (distances / (flowspeed * 3.6) + 0.5).astype(int)

        # Aggregation of all basin flows (m^3/hr):
        # Calculating the  timedelayed inflow at each plant
        # The the offset is done by cyclic shift in hours
        for b in ppl.upstream:
            inflow_plant += basin_inflows.sel(hid=b).roll(time=nhours.at[b]) # MODFIED += doesn't work...

    return inflow.to_pandas().T

def calculate_ror_power(sites, inflow_series):
    '''
    This calculate the power avaialbiliy series based on a potential energy approach
    using the inflow series, power capacity, and annual energy production for the facility.
    The annual energy production is used to normalize the series. In this implementation Newton's
    Method is used to find the constant which satisfies the following equation:
    E_tot = sum_t(min(const*inflow_t, capacity_rating))
    This function finds the "const" the make the equality above true.
    '''
    n_itr = 100 # gives 100 iterations to converge
    power_series = inflow_series.copy() # Store power availability series
    for ind,row in sites.iterrows():
        # initialize value for first points of evaluation
        const = row['annual_avg_energy'] * 1000 / inflow_series[ind].sum()
        for i in range(n_itr):
            # Calculate power series based on const
            try:
                temp_series = inflow_series[ind].apply(lambda x: min(x*const,row['capacity']))
            except:
                print('Error in truncated power series in calculate_ror_power')
                exit(1)
            # Calculate function value (this is what we are trying to get to zero)
            f = -row['annual_avg_energy'] * 1000 + temp_series.sum()
            if abs(f) <= 1e-1:
                power_series[ind] = temp_series
                break
            # Calculate derivative
            try:
                dfdc = inflow_series[ind].apply(lambda x: (x*const < row['capacity'])*x)
            except:
                print('Error in computing dfdc in calculate_ror_power')
                exit(1)
            if dfdc.sum() == 0:
                const = const / 2 # make const smaller as likely initilized too large 
            else:
                # Perform newton method update
                const = const - f / dfdc.sum()
            # conditon check
        if i == (n_itr-1): # Convergence
            print('Newton Method did not converge for RoR power series.')
            print(f'Final function value is {f}')
            exit(2)

    return power_series


def create_ror_power(sites_prep, basin_data, cutout, cfg):
    '''
    This function is used to calculate RoR power based on the potential energy approach.
    '''
    # Reduce assets to only asset_id since asset_id is 1-to-1 with inflows
    subset =["asset_id", "lat", "lon"]
    sum_list = ["capacity", "annual_avg_energy"]
    sites = sites_prep.groupby(by="asset_id", group_keys=False).apply(lambda x:
                                                                     merge_assets(x, subset, sum_list))

    # (i) Calculated the inflows for each site
    basins = prepare_basins(sites, basin_data)
    basin_inflows = calculate_basin_inflows(basins, cutout, height=bool(cfg['ror_inflows']['height']))
    site_inflows = calculate_plant_inflows(basin_inflows, basins, flowspeed=cfg['ror_inflows']['flowspeed'])

    # (ii) Compute the power availability series for each site based on capacity, inflow series, and
    # annual energy production
    power_series = calculate_ror_power(sites, site_inflows)
    power_series.to_csv(cfg['ror_inflows']['ror_outfile'])


def calibrate_reservoir_inflow(site_inflows, fpath, method="mean_inflow_calibrate"):
    '''
    fpath: Path to location of the inflow tables with statistics of historical inflows
           for BC reservoirs.
    '''
    # loop over reservoirs to normalize\
    # Final flow is in flow per hour
    
    inflows_adjusted = pd.DataFrame(index=pd.to_datetime([]))
    for rid in site_inflows.columns: # loop over reservoirs

        if method == "mean_inflow_calibrate":
            inflows_adjusted[rid] = mean_inflow_calibration(site_inflows[rid], rid, fpath)

        elif method == "mean_inflow_constant":
            inflows_adjusted[rid] = mean_inflow_constant(site_inflows[rid], rid, fpath)

        else:
            print(f'The method {method} for calculating reservoir inflows has not been created')

    return inflows_adjusted

def mean_inflow_calibration(site_inflows, rid, fpath):
    '''
    This function is designed to take in a unnormalized inflow series and normalize it.
    The normalize is done for each month based on the mean monthly inflow in (cms).
    This function is applied to a single month segment at a time.
    '''
    month_2_num = {"January":1, "February":2, "March":3, "April":4, "May":5,
                        "June":6, "July":7, "August":8, "September":9,
                        "October":10, "November":11, "December":12}

    temp_series = pd.Series(dtype='float64')
    read_path = fpath + f"/{rid}" + ".csv" # Format should be hardcoded in a configuration file

    if not os.path.exists(read_path):
        print(f'There are no flow statistics for {rid}')
        exit(123)

    inflow_stats = pd.read_csv(read_path)

    for idx,row in inflow_stats.iterrows(): # loop over monthly stats 
        # prep mask
        mask = (site_inflows.index.month == month_2_num[row['Month']])

        # calculations
        inflow_seg = site_inflows.loc[mask]
        num_hours = inflow_seg.shape[0]
        const = inflow_seg.sum() / (3600 * num_hours * row['Mean Monthly Inflow'])
        norm_inflow_seg = inflow_seg / const

        # store results (month-by-month)
        temp_series = pd.concat([temp_series, norm_inflow_seg])
        
    return temp_series.sort_index()


def mean_inflow_constant(site_inflows, rid, fpath):
    '''
    This function is designed to return a constant monthly inflow equal to
    that reservoirs historical mean inflow for each month.
    '''
    month_2_num = {"January":1, "February":2, "March":3, "April":4, "May":5,
                        "June":6, "July":7, "August":8, "September":9,
                        "October":10, "November":11, "December":12}

    temp_series = pd.Series(dtype='float64')
    read_path = fpath + f"/{rid}" + ".csv"

    if not os.path.exists(read_path):
        print(f'There are no flow statistics for {rid}')
        exit(123)

    inflow_stats = pd.read_csv(read_path)

    # loop over df_inflow_stats for this particular reservoir (month by month)
    for idx,row in inflow_stats.iterrows():
        mask = (site_inflows.index.month == month_2_num[row['Month']]) 
        # calculations
        inflow_seg = site_inflows.loc[mask]

        norm_inflow_seg = pd.Series([row['Mean Monthly Inflow'] * 3600] * inflow_seg.shape[0],
                                dtype='float64', index=inflow_seg.index)


        temp_series = pd.concat([temp_series, norm_inflow_seg])

    return temp_series.sort_index()


def create_cascade_inflow(reservoir_sites, basin_data, cutout, hydro_sites, cfg, method = "mean_inflow_calibrate"):
    '''
    This function creates the file containing the inflows for each reservoir apart of the cascades.
    There are 2 approaches used for calculating the inflows based on following conditions
    1) If a reservoir has a water use plan with monthly statistics. Then a normalization approach is used
        normalizing/calibrating with the monthly statistic.
    2) If a reservoir does not have a water use plan with monthly statistics AND it is a downstream reservoir
        then it is imputed with an inflow time_series of 0.
    reservoir_sites: Load all reservoir assets.
    basin_data: GeoDataFrame of all the hydrobasin data.
    cutout: Cutout used by Atlite for merging inflows into basins.
    hydro_sites: GeoDataFrame of hydro generation sites.
    cfg: Dictionary with all configuration options.
    method: Method used to normalize the data.
    '''

    basins = prepare_cascade_basins(reservoir_sites, basin_data, hydro_sites)

    basin_inflows = calculate_basin_inflows(basins, cutout,
                                                height=bool(cfg['reservoir_inflows']['height']))

    site_inflows = calculate_plant_inflows(basin_inflows, basins,
                                                flowspeed=cfg['reservoir_inflows']['flowspeed'])
    
    # File path + name for reading in inflow tables
    fpath = cfg['bc_hydro']['inflow_tables']

    # 1) normalize inflow time series for selected reservoirs
    reservoirs = hydro_sites[hydro_sites['hydro_type'] == "reservoir"]['upper_reservoir_id'].unique().tolist() # TRY Unique it
    final_inflows = calibrate_reservoir_inflow(site_inflows[reservoirs], fpath, method)

    # 2) Currently impute with for for reservoirs downstream with no WUP statistics.
    reservoirs_impute = hydro_sites[hydro_sites['hydro_type'] == "reservoir-impute"]['upper_reservoir_id'].tolist()
    for rid in reservoirs_impute:
        final_inflows[rid] = pd.Series([0]*final_inflows.shape[0],
                                     dtype='float64', index=final_inflows.index)
        

    # power_series = calculate_ror_power(sites, site_inflows)
    final_inflows.to_csv(cfg['reservoir_inflows']['output'])


def check_wup_exists(rid, fpath):
    '''
    Checks if water use plan statistics (monthly) exist for a particular reservoir
    return: boolean (true if it exists)
    '''
    # Format should be in a configuration file or something similar eventually
    read_path = fpath + f"/{rid}" + ".csv" 

    return os.path.exists(read_path)

def check_head_reservoir(rid, cascade_sites):
    '''
    This function deteremines if a particular rid is a downstream_reservoir_id for any other sites.
    If it is not then that means the rid is an head upstream reservoir (i.e. not inflows caused by other reservoirs)

    rid: upper_reservoir_id.
    cascade_sites: dataframe of all cascade sites. Specifically, those which do not have WUP inflow statistics.
    return: boolean. Indicating whether the particular reservoir is a head reservoir (True) otherwise (False).
    '''
    return not (cascade_sites['lower_reservoir_id'] == rid).any()

# Code below this point is used for formatting of the hydroelectric data


if __name__ == '__main__':
    print('This module contains all preprocessing functions used for creating all features needed to model Hydroelectric facilities in PyPSA_BC')
