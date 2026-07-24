import pandas as pd
import plotly.express as px

def create_generation_plots(df:pd.DataFrame):
    """
    Create an interactive area plot for generator power output over time.

    Parameters:
    df (pd.DataFrame): DataFrame containing the power output data with time index.

    Returns:
    fig: Plotly figure object representing the area plot.
    """
    # Create an interactive line plot
    fig = px.area(df, template="plotly_white", title="Simulated generation (2040)", labels={"value": "Power", "index": "Time"})
    fig.update_layout(xaxis_title="Time", yaxis_title="MW")
    return fig