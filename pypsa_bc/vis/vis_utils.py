import pandas as pd
from pathlib import Path
from typing import Optional
import plotly.express as px
from pypsa_bc import utils
import plotly.graph_objects as go

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



def create_generation_plots(df:pd.DataFrame):
    """
    Create an interactive area plot for generator power output over time.

    Parameters:
    df (pd.DataFrame): DataFrame containing the power output data with time index.

    Returns:
    fig: Plotly figure object representing the area plot.
    """
    # Create an interactive line plot

    fig = px.area(df, template="plotly_white", title="Simulated generation", labels={"value": "Power", "index": "Time"})
    fig.update_layout(xaxis_title="Time", yaxis_title="MW")
    return fig


def plot_hourly_profile(data: pd.DataFrame):
    title = f'Normalized Load Profile {data.index[0].year}'
    peak_load_time = data['LOAD'].idxmax()
    peak_load_value = data['LOAD_profile_norm_peak2hr'].max()

    # Function to resample data
    def resample_data(freq):
        return data.resample(freq).mean()

    # Create the figure
    fig = go.Figure()

    # Add hourly data as the main trace
    fig.add_trace(go.Scatter(
        x=data.index,
        y=data['LOAD_profile_norm_peak2hr'],
        mode='lines',
        name='Hourly',
        line=dict(color='skyblue')
    ))

    # Add resampled data for different frequencies
    freq_options = ['D', 'W', 'ME']  # Daily, Weekly, Monthly
    freq_options_name = {
        'D': 'Daily',
        'W': 'Weekly',
        'ME': 'Monthly'
    }

    for freq in freq_options:
        resampled_data = resample_data(freq)
        fig.add_trace(go.Scatter(
            x=resampled_data.index,
            y=resampled_data['LOAD_profile_norm_peak2hr'],
            mode='lines',
            name=freq_options_name[freq],
            line=dict(dash='dot')  # Dotted line style
        ))

    # Add peak load as a red circle
    fig.add_trace(go.Scatter(
        x=[peak_load_time],
        y=[peak_load_value],
        mode='markers',
        marker=dict(color='red', size=12, symbol='circle'),
        name='Peak Load'
    ))

    # Customize layout
    fig.update_layout(
        template='plotly_white',
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18}
        },
        xaxis_title=None,
        yaxis_title="Normalized Load Profile",
        xaxis=dict(showgrid=False, ticks="outside"),
        yaxis=dict(
            showgrid=True,
            gridcolor="lightgrey",
            tickformat='.2f'  # Format y-axis as decimals
        ),
        font=dict(size=14),
        margin=dict(l=40, r=40, t=50, b=40),
        showlegend=True
    )

    return fig