import sys
from pathlib import Path
import pandas as pd
import geopandas as gpd
import shapely as sp
from bc_power import disaggregation_utils

#Get dataframe of centers depending on what was requested (geograpic center or population center)
def get_center(gp, bc, muni):
    regions = bc['NAME_2'].rename('REGION')

    if gp == 'geo':
        geocenter = bc['geometry'].apply(lambda x: x.centroid.coords[0]).reset_index()
        geocenter = geocenter.drop(columns='index').rename(columns={'geometry': 'geocenter'})

        geo_lon = geocenter.apply(lambda x: x[0][0], axis=1)
        geo_lat = geocenter.apply(lambda x: x[0][1], axis=1)

        centers = pd.DataFrame(data={'REGION': regions, 'LONGITUDE': geo_lon, 'LATITUDE': geo_lat}).set_index('REGION')

    elif gp == 'pop':
        #Empty dataframe to be filled in the loop below
        centers = gpd.GeoDataFrame(columns=['LONGITUDE', 'LATITUDE'])

        #Loop through each region
        for region in muni['region'].unique():
            #Get a subset of the dataframe containing just the municipalities for one region
            subset = muni.loc[muni['region'] == region].copy()

            #Total population to use for calculating population center
            total_population = subset['population'].sum()

            #Weighted average of the cartesian coordinates for municipalities
            subset['cartesian'] = subset.apply(lambda x: disaggregation_utils.to_cartesian(x['latitude'], x['longitude']), axis=1) #Get cartesian coordinates

            avg_x = subset.apply(lambda x: x['cartesian'][0] * x['population'], axis=1).sum() / total_population #x = [0]
            avg_y = subset.apply(lambda x: x['cartesian'][1] * x['population'], axis=1).sum() / total_population #y = [1]
            avg_z = subset.apply(lambda x: x['cartesian'][2] * x['population'], axis=1).sum() / total_population #z = [2]

            lat, lon = disaggregation_utils.to_lat_lon(avg_x, avg_y, avg_z)
            
            #Add this point to pop_center
            add = gpd.GeoDataFrame(data={'LONGITUDE': [lon], 'LATITUDE': [lat]})
            centers = pd.concat([centers, add])

        #Fixing index
        centers = centers.reset_index()
        centers = centers.drop(columns='index').set_index(regions)
    
    else:
        centers = None

    return centers

#Function to assign a substation to a region in BC based on its longitude and latitude
def get_region(row, bc):
    #Define point p as the (longitude, latitude) for a substation
    p = sp.geometry.Point(row.longitude, row.latitude)

    #Flags to check if a substation is within a region's polygon
    region = bc[p.within(bc.geometry)]
    
    #Return None if a station cannot be pinned down to any particular region
    if len(region) == 0:
        return None
    
    #Return the region that the substation is located at
    return region['NAME_2'].values[0]

#Disaggregate to substation level here
def disaggregate(bc, res, csmi, substations, centers):
    #Assign all of the substations to a region within BC
    bc_stations = substations.copy()
    bc_stations['region'] = substations.apply(get_region, bc=bc, axis=1)
    bc_stations = bc_stations.dropna() #Need to drop None values just in case some substations weren't able to be mapped to a region

    for region in bc['NAME_2'].unique():
        #Subset of bc_stations by region
        subset = bc_stations.loc[bc_stations['region'] == region].copy()

        #The region's center 'c' as a point [x, y, z]
        c = disaggregation_utils.to_cartesian(centers.loc[region]['LATITUDE'], centers.loc[region]['LONGITUDE'])

        #Get cartesian coordinates for all of the substations within the region
        subset['cartesian'] = subset.apply(lambda x: disaggregation_utils.to_cartesian(x['latitude'], x['longitude']), axis=1)

        #Need a separate subset just for residential, as it does not receive power from ISS nodes
        subset_res = subset.loc[subset['node_type'].isin(['DSS', 'DFS', 'TSS'])].copy()

        #This check is needed because it is possible for a region to only have ISS nodes
        if len(subset_res) > 0:
            #Proportion calculations for residential
            subset_res['proportion'] = subset_res.apply(lambda x: disaggregation_utils.dist_3d(c, x['cartesian']), axis=1)
            res_total = subset_res['proportion'].sum()
            subset_res['proportion'] = subset_res['proportion'].apply(lambda x: x / res_total)

            print(res[region])


    return 0, 2



#Does some input verification before generating the regional hourly loads
'''

    USAGE:

    python dissagregate_load.py <ARG1> <ARG2> <ARG3> <ARG4>

    ARG1 = The CEEI spreadsheet file
    ARG2 = The folder which contains the hourly load data from BC Hydro. This folder should contain files named BalancingAuthorityLoad20XX.xls
    ARG3 = The year to use for disaggregation
    ARG4 = The folder to write the outputs to. This script outputs two CSV files, one for residential load, one for industrial load

    EXAMPLE:

    python disaggregate_load.py CEEI_2020.xlsx .\load\ 2015 .\hourly_load_region\
    
'''
def main():
    try:
        #Try reading the arguments passed in the terminal
        gadm_path = Path(sys.argv[1])
        hourly_res_path = Path(sys.argv[2] + '/hourly_res_' + sys.argv[3] + '.csv')
        hourly_csmi_path = Path(sys.argv[2] + '/hourly_csmi_' + sys.argv[3] + '.csv')
        year = int(sys.argv[3])
        substations_path = Path(sys.argv[4])
        center = sys.argv[5] # 'geo' or 'pop'
        municipal_path = Path(sys.argv[6])
        output_path_res = Path(sys.argv[7] + '/substation_hourly_res_' + sys.argv[3] + '.csv')
        output_path_csmi = Path(sys.argv[7] + '/substation_hourly_csmi_' + sys.argv[3] + '.csv')

    except Exception as e:
        #Less than 7 arguments, return error code 1
        print('There are inputs missing')
        print(e)
        return 1
    
    else:
        #Correct number of arguments
        try:
            #Try loading in the data
            gadm = gpd.GeoDataFrame.from_file(gadm_path)
            bc = gadm.loc[gadm.NAME_1 == 'BritishColumbia'].reset_index()
            bc = bc.drop(columns='index') #Filtered down to just the subdivisions of BC

            #Hourly load data for the subdivisions of BC
            hourly_res = pd.read_csv(hourly_res_path, parse_dates=True, index_col='TIME')
            hourly_csmi = pd.read_csv(hourly_csmi_path, parse_dates=True, index_col='TIME')

            #Substations of BC
            substations = pd.read_csv(substations_path)
            bc_stations = substations.loc[(substations.province == 'BC') & (substations.node_type.isin(['DSS', 'ISS', 'DFS', 'TSS']))].copy().drop(columns=['owner', 'province', 'operating_region', 'sources', 'notes'])

            bc_stations['region_polygon'] = bc_stations.apply(get_region, bc=bc, axis=1)
            bc_stations.to_csv('../../regions_and_substations.csv', index=False)

            #Municipality populations of BC
            municipal = pd.read_csv(municipal_path, index_col='index')

            #Centers either geographic or population centers
            centers = get_center(center, bc, municipal)

            if centers is None:
                raise ValueError("Check if you entered 'geo' or 'pop' in the arguments")


        except Exception as e:
            #An input may be spelled incorrectly, return error code 2
            print('One or more inputs are in the wrong format')
            print(e)
            return 2
        else:
            #All is good, start distributing to substations

            substations_res, substations_csmi = disaggregate(bc, hourly_res, hourly_csmi, bc_stations, centers)

            #Write to files to the output folder path
            #substations_res.to_csv(output_path_res); substations_csmi.to_csv(output_path_csmi)

            print(substations_res)
            print(substations_csmi)
            
            #Return code 0 is for when everything runs without a problem
            return 0
        

if __name__ == '__main__':
    main()