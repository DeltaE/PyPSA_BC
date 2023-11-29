import atlite
import pandas as pd
import geopandas as gpd
import xarray as xr
from bc_power import wind

'''
GENERATING POWER TIME SERIES HERE
'''

#Function to calculate the generation profiles from the cutout, for a given farm
#Used in generate_wind_ts(), generate_solar_ts()
#cutout = Cutout object
#cluster = Generator cluster, these are rows from X_assets.csv
def calculate_MW(cutout, cluster, type):
    '''
    This function calculates a timeseries of power production in MW.
    cutout: Atlite cutout object
    cluster: A specific VRE asset (wind farm or solar-pv farm), row of a dataframe.
    type: Type of VRE, "wind" or "pv" are the only solutions supported right now.
    '''
    if type == 'wind':
        cap_factors = cutout.wind(turbine=wind.get_config(cluster['config_oedb'], cluster['Hub height (m)']), capacity_factor=True)
    elif type == 'solar': # PV Panel
        cap_factors = cutout.pv(panel='CdTe', orientation='latitude_optimal', capacity_factor=True)
    else:
        print(f"error in the calculate_MW function type:{type} not implemented yet!")
        exit(0)
    # This code calculates the cells/grids for where generation exists
    # Cells of the cutout
    cells = cutout.grid 

    #site: location of the turbines with latitudes and longitudes specified
    sites = gpd.GeoDataFrame([[cluster['location'], cluster['longitude'], cluster['latitude'], cluster['Install capacity']]],
                            columns=['name', 'lon', 'lat', 'capacity']
                            ).set_index('name')

    # Finds cutout cells nearest to (lat,lon) of each site 
    # (x,y) of cells/grids are the center points by default and correspong to (lat,lon) of the cells.
    nearest = cutout.data.sel(
        {'x': sites.lon.values, 'y': sites.lat.values}, 'nearest').coords

    # Add new columns for the grid (x,y) where each generator falls.
    sites['x'] = nearest.get('x').values
    sites['y'] = nearest.get('y').values

    # Performs an inner join on the (x,y) values
    cells_generation = sites.merge(
        cells, how='inner').rename(pd.Series(sites.index))

    # Layout indicates the installed capacity and grid location of installed capacity
    layout = xr.DataArray(cells_generation.set_index(['y', 'x']).capacity.unstack())\
                        .reindex_like(cap_factors).rename('Installed Capacity [MW]')

    # Finally generate the power produced at these sites
    # Returns the power generation of the generators for each location
    # The data is returned as an xarray with shape of (timesteps, generators)
    if type == 'wind':
        power_generation = cutout.wind(turbine=wind.get_config(cluster['config_oedb'], cluster['Hub height (m)']),
                                        layout = layout,
                                        shapes = cells_generation.geometry,
                                        per_unit = False) # currently using per-unit
    elif type == 'solar':
        power_generation = cutout.pv(panel='CdTe', 
                                        orientation='latitude_optimal',
                                        layout = layout,
                                        shapes = cells_generation.geometry,
                                        per_unit = False) # currently using per-unit
    else:
        print(f"error in the calculate_MW function type:{type} not implemented yet!")
        exit(0)

    return power_generation.to_pandas()

'''
CALIBRATING TIME SERIES TO CODERS AAG HERE
'''

#Function to calculate some value 'a' such that sum(a*w) = P_rated
#Used in calibrate_generation()
#w_vec = Vector of generation values
#aag = Target AAG to calibrate to
#p = P_rated
#a_init = Initial guess for the value of 'a'
def get_A(w_vec, aag, p, a_init):
    #Initialize values for the loop
    iter = 0
    a = a_init

    while iter < 100: #Stop after 100 iterations if it does not work out for some reason (this should never happen in practice)
        #dfda:
        #-w, when 0 <= a*w <= P_rated
        #0, everywhere else
        dfda = w_vec.apply(lambda w: float(a*w >= 0) * float(a*w <= p) * w * -1).sum()

        #The value for f = AAG - sum of t_{i} where t_{i} = min(a*w_i, P_rated)
        f = aag - w_vec.apply(lambda w: min(a*w, p)).sum()

        #Get new value for 'a', a_{i+1} = a_{i} - (f / dfda)
        a_new = a - (f / dfda)

        #Calculate tolerance percentage
        tol = abs(100 * ((a_new - a) / a))

        #What to do based on the value of the tolerance
        if tol < 0.01:
            break #Exit loop and return aNew if tolerance is small enough
        else:
            a = a_new #Loop again using the newly calculated 'a' value

        iter += 1
    
    return a_new

#Function to calibrate the wind generation data by using the CODERS AAG data
#Used in main()
#wind_ts = Wind time series data frame from generate_wind_ts()
#wind_assets = Wind assets data frame from wind_assets.py, this has the CODERS AAG data ready
def calibrate_generation(ts, assets):
    #Get date breakpoints for each year in the wind time series
    years = ts.groupby(ts.index.year).count()[ts.columns[0]]
    year_index = [0]

    #Build year_index, which will be used to calibrate wind generations per year
    for year in years.index:
        year_index.append(year_index[-1] + years.loc[year])
    
    #Empty DataFrame whose columns will be filled
    calibrated = pd.DataFrame(data={'time': ts.index})

    #Generate calibrated wind generation time series
    for component in assets['asset_id'].unique(): #NOTE: component_id -> asset_id (2023-11-15)
        #Get necessary parameters for calculating the value 'a' for calibration
        gen = ts[component].reset_index()
        aag = assets.loc[assets['asset_id'] == component]['CODERS AAG (GWh/y)'].values[0] * 1000 #GWh --> MWh
        p = assets.loc[assets['asset_id'] == component]['Install capacity'].sum()
        
        #Empty series to hold the generation for a component_id
        column_to_add = pd.Series(name=component, dtype='float64')

        #Calculate calibrated wind generations by year
        for i in range(len(year_index) - 1):
            subset = gen.loc[(gen.index >= year_index[i]) & (gen.index < year_index[i+1])].copy()

            #Obtain a value of a for the year
            # NOTE: This function to calibrate should be standardized to be same as run-of-river
            a = get_A(subset[component], aag, p, aag / subset[component].sum())

            #Use the a value to normalize the wind generation in that year
            subset[component] = subset[component].apply(lambda w: min(a*w, p))

            column_to_add = pd.concat([column_to_add, subset[component]], ignore_index=True)
        
        #Add the column to the final data drame
        calibrated = pd.concat([calibrated, column_to_add], axis=1)
    

    #Final calibrated generation data frame to return, set index to time to match the uncalibrated data frame's format
    return calibrated.set_index('time')