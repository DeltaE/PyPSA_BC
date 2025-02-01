import pandas as pd
from pathlib import Path
from typing import Optional
import plotly.express as px
from pypsa_bc import utils

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
                         type='area',
                         plot_title:str=None,
                         yaxis_title:str=None,
                         xaxis_title:str=None,
                         legend_title_text:str=None,
                         save_to:Optional[Path]=None,
                         show:Optional[bool]=True):
    if type=='area':
        fig = px.area(data, title="Title" if plot_title is None else plot_title)
    if type=='line':
        fig = px.line(data, title="Title" if plot_title is None else plot_title)
    
    fig.update_layout(xaxis_title='Time' if xaxis_title is None else xaxis_title , 
                      yaxis_title='Variable' if yaxis_title is None else yaxis_title,
                      legend_title_text="" if legend_title_text is None else legend_title_text)
  
    if save_to is not None:
        save_to = Path(save_to)
        save_to.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(save_to)
        utils.print_update(level=2,message=f'Plot save to :{save_to}')
    if show:
        fig.show()