import pandas as pd
from pathlib import Path
from typing import Optional
import plotly.express as px
from pypsa_bc import utils
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import contextily as ctx
import geopandas as gpd
import matplotlib as mpl
import pypsa
from pypsa_bc import calc

"""
def load_and_process_data(file_path, 
                          technologies):
    df = pd.read_csv(file_path)
    year_col = 'YEAR'
    technology_col = 'TECHNOLOGY'
    if 'YEAR' not in df.columns or 'TECHNOLOGY' not in df.columns:
        year_col = 'y'
        technology_col = 't'
    yearly_summed_values = {tech: [] for tech in technologies}
    all_years = sorted(set(df[year_col]))
    grouped_data = {}
    for tech in technologies:
        filtered_df = df[df[technology_col].str.startswith(tech)]
        last_column = filtered_df.iloc[:, -1]
        grouped = last_column.groupby(filtered_df[year_col]).sum()
        grouped_data[tech] = grouped.reindex(all_years, fill_value=0)
        yearly_summed_values[tech] = [grouped.get(year, 0) for year in all_years]
    result_df = pd.DataFrame(grouped_data)

    return all_years, result_df

""" 

def load_and_process_data(file_path, 
                          technologies):
    df = pd.read_csv(file_path)
    year_col = 'YEAR'
    technology_col = 'TECHNOLOGY'
    if 'YEAR' not in df.columns or 'TECHNOLOGY' not in df.columns:
        year_col = 'y'
        technology_col = 't'
    yearly_summed_values = {tech: [] for tech in technologies}
    all_years = sorted(set(df[year_col]))
    grouped_data = {}
    for tech in technologies:
        filtered_df = df[df[technology_col].str.startswith(tech)]
        last_column = filtered_df.iloc[:, -1]
        grouped = last_column.groupby(filtered_df[year_col]).sum()
        grouped_data[tech] = grouped.reindex(all_years) #, fill_value=0
        yearly_summed_values[tech] = [grouped.get(year, 0) for year in all_years]
    result_df = pd.DataFrame(grouped_data)

    return all_years, result_df

def visualize_timeseries(data:pd.DataFrame,
                         year:int,
                         type='area',
                         plot_title:str=None,
                         yaxis_title:str=None,
                         xaxis_title:str=None,
                         legend_title_text:str=None,
                         save_to:Optional[Path]=None,
                         show:Optional[bool]=False):
    data.index = data.index.map(lambda x: x.replace(year=year))
    if type=='area':
        fig = px.area(data, title="Title" if plot_title is None else plot_title)
    if type=='line':
        fig = px.line(data, title="Title" if plot_title is None else plot_title)
    
    fig.update_layout(xaxis_title=None if xaxis_title is None else xaxis_title , 
                      yaxis_title='Variable' if yaxis_title is None else yaxis_title,
                      legend_title_text="" if legend_title_text is None else legend_title_text)
  
    if save_to is not None:
        save_to = Path(save_to)
        save_to.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(save_to)
        utils.print_update(level=2,message=f'Plot save to :{save_to}')
        
    if show:
        fig.show()
        
    return fig

def plot_inter_region_link_usage(year:int,
                                 pypsa_network:pypsa.Network,
                                 regional_boundaries_GADM_L2:gpd.GeoDataFrame,
                                 usage_statistics='usage_avg',
                                 plot_dpi:int=300,
                                 line_cmap:str='Reds',
                                 plot_title=None,
                                 show:bool=False,
                                 plot_save_to:str|Path=None):
    
    usage_type='Average' if 'avg' in usage_statistics else "Max"
    title=f'Inter Region Link {usage_type} Usage [Simulated for {year}]' if plot_title is None else plot_title

    inter_region_line_with_usage=calc.get_line_usage(pypsa_network,regional_boundaries_GADM_L2)

    if 'Region' not in regional_boundaries_GADM_L2.columns:
        regional_boundaries_GADM_L2 = regional_boundaries_GADM_L2.reset_index(inplace=True)
    
    boundary = regional_boundaries_GADM_L2
    
    # Assign a number to each region
    boundary['Region_Number'] = range(1, len(boundary) + 1)

    # Transform GeoDataFrames to match basemap CRS (EPSG:3857)
    boundary = boundary.to_crs(epsg=3857)
    inter_region_lines = inter_region_line_with_usage.to_crs(epsg=3857)
    
    # Scale s_nom to line widths between 1 and 10
    min_width, max_width = 2, 10
    s_nom = inter_region_line_with_usage['s_nom']
    scaled_width = (s_nom - s_nom.min()) / (s_nom.max() - s_nom.min())  # Normalize to 0-1
    inter_region_line_with_usage['line_width'] = min_width + scaled_width * (max_width - min_width)
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot boundary and inter-region lines
    boundary.plot(ax=ax, color='None', edgecolor='k', linewidth=0.3, alpha=0.5)

    # Set the bounds for the colorbar (can adjust these based on data range)
    vmin, vmax = 0, 1

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    inter_region_lines.plot(column=usage_statistics, 
                            cmap=line_cmap,
                            linewidth=inter_region_line_with_usage['line_width'],
                            alpha=1,
                            legend=False,  # We'll add a custom legend for the colorbar
                            ax=ax,
                            norm=norm)

    # Annotate region numbers at centroids
    for idx, row in boundary.iterrows():
        centroid = row.geometry.centroid
        ax.annotate(f"{row['Region_Number']}", 
                    xy=(centroid.x, centroid.y), 
                    horizontalalignment='center', fontsize=8, color='black',
                    bbox=dict(facecolor='none', edgecolor='none'))

    # Add basemap
    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)

    # Turn off the grid and axis
    ax.axis('off')

    # Create a ScalarMappable for the custom colorbar
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=line_cmap, norm=norm)
    sm.set_array([])  # The array is empty since we're just setting the colorbar

    # Add colorbar with custom settings
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', shrink=0.8, pad=0.001, anchor=(.5, .8))
    cbar.set_label(f'{usage_type} Usage', rotation=0, labelpad=5,fontsize=12)
    cbar.outline.set_visible(False)  # Remove colorbar border
    cbar.ax.patch.set_alpha(0.7)  # Add transparency to colorbar

    # Create custom legend for region number mapping
    handles = [
        Line2D([0], [0], marker=None, color='w', markerfacecolor=None, markersize=None,
               label=f"{row['Region_Number']}: {row['Region']}",
               ) 
        for idx, row in boundary.iterrows()
    ]

    # Place the legend outside the plot
    plt.legend(handles=handles, loc='upper left', bbox_to_anchor=(1, .98), frameon=False)

    # Add title
    plt.title(title)

    # Save the plot
    plt.tight_layout()
    if plot_save_to:
        plot_save_to = Path(plot_save_to)
        plot_save_to.parent.mkdir(exist_ok=True, parents=True)
        plt.savefig(plot_save_to, dpi=plot_dpi)
    plt.close(fig)
    return fig