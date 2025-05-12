import pandas as pd
import plotly.express as px

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def create_generation_plots(df:pd.DataFrame,
                            year:int):
    """
    Create an interactive area plot for generator power output over time.

    Parameters:
    df (pd.DataFrame): DataFrame containing the power output data with time index.

    Returns:
    fig: Plotly figure object representing the area plot.
    """
    # Create an interactive line plot
    df.index = df.index.map(lambda x: x.replace(year=year))
    fig = px.area(df, template="plotly_white", title="Simulated generation (2040)", labels={"value": "Power", "index": "Time"})
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
