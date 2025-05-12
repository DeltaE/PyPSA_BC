from shapely.geometry import LineString
import pypsa
import geopandas as gpd
import pandas as pd

def get_line_usage(solved_network:pypsa.Network,
                   regional_boundaries_GADM_L2:gpd.GeoDataFrame):
    
    n=solved_network
    
    regional_boundaries_GADM_L2['Region'] = regional_boundaries_GADM_L2['Region'].str.replace(' ', '')
    boundary=regional_boundaries_GADM_L2.set_index('Region')
    
    # Creates Boundary geoms mapping to region-bus
    n.lines['bus0_geom'] = n.lines.bus0.apply(lambda x: boundary.loc[x].geometry if x in boundary.index else None)
    n.lines['bus1_geom'] = n.lines.bus1.apply(lambda x: boundary.loc[x].geometry if x in boundary.index else None)

    
    # Pre-allocate the list for geometries
    geometries:list = []

    # Vectorized check for non-null 'bus0_geom' and 'bus1_geom'
    valid_rows = n.lines['bus0_geom'].notnull() & n.lines['bus1_geom'].notnull()

    # Only process the valid rows
    for row in n.lines[valid_rows].itertuples():
        line = LineString([row.bus0_geom.centroid, row.bus1_geom.centroid])
        geometries.append(line)

    # For invalid rows, append None
    geometries.extend([None] * (len(n.lines) - len(geometries)))

    # Create the GeoDataFrame with the valid geometries
    inter_region_lines:gpd.GeoDataFrame = gpd.GeoDataFrame(n.lines, geometry=geometries, crs='4326')

    # Drop rows with None geometries (invalid lines)
    inter_region_lines = inter_region_lines[inter_region_lines.geometry.notnull()]

    link_usage:list= []
    
    for link in n.lines_t.p0.columns:
        max_value = n.lines_t.p0[link].apply(lambda x: x / n.lines.loc[link, 's_nom']).max()
        avg_value = n.lines_t.p0[link].apply(lambda x: x / n.lines.loc[link, 's_nom']).mean()
        # min_value = n.lines_t.p0[link].apply(lambda x: x / n.lines.loc[link, 's_nom']).min()
        # std_dev = n.lines_t.p0[link].apply(lambda x: x / n.lines.loc[link, 's_nom']).std()
        link_usage.append({
            'Link': link, 
            'Max': abs(max_value), 
            'Average': abs(avg_value), 
            # 'Min': abs(min_value), 
            # 'StdDev': abs(std_dev)
        })

    link_usage_df = pd.DataFrame(link_usage)
    link_usage_df=link_usage_df.sort_values(by='Max',ascending=False)
    link_usage_df.set_index('Link',inplace=True)
    inter_region_lines['usage_max'] = inter_region_lines.index.map(link_usage_df['Max'])
    inter_region_lines['usage_avg'] = inter_region_lines.index.map(link_usage_df['Average'])
    inter_region_lines=inter_region_lines.sort_values(by='usage_max',ascending=False)
    inter_region_lines
    return inter_region_lines