from dash import Dash, html, dcc, Output, Input
import pandas as pd
from bc_combined_modelling.vis import vis_utils
from bc_combined_modelling import utils
import plotly.express as px
from bc_combined_modelling.vis.dashboard.dash_components import header, scenario_and_tabs
from pathlib import Path
import plotly.graph_objects as go
from bc_combined_modelling.vis.dashboard.plot_demand import create_demand_plots,create_demand_plot_simplified
from bc_combined_modelling.vis.dashboard.plot_pypsa_gen import create_generation_plots
# Load configurations
dash_configs: dict = utils.load_config('config/dashboard.yaml')

plot_files: dict = dash_configs['plot_files']
demand_data = dash_configs['demand_data']
units_mapping = dash_configs['units_mapping']
filenames_mapping = dash_configs['filenames_mapping']
selected_technologies = dash_configs['technologies']
custom_colors = dash_configs['custom_colors']
legend_labels = dash_configs['legend_labels']
result_type = dash_configs['result_type']
info_SCENARIOS = dash_configs['info_SCENARIOS']

# Define folders/sub-folders to access data
scenario_results_directory = Path('results/clews')
SCENARIOS = sorted([item.name for item in scenario_results_directory.iterdir() if item.is_dir()])
# Add 'Resources' to plot_files and SCENARIOS if required
plot_files['resources'] = []  # No specific files needed for this tab
plot_files['dispatch'] = []  # No specific files needed for this tab


timeslices = 8
results_subfolder = f'{timeslices}ts_csvs_gurobi'
""" 
folder_names = []
for scenario in SCENARIOS:
    # Define the path to the directory
    directory_path = f'results/clews/{scenario}'  # Replace with your actual path
    if os.path.exists(directory_path):
        # List all subdirectories that end with 'ts'
        folder_names.extend([folder for folder in os.listdir(directory_path) 
                             if os.path.isdir(os.path.join(directory_path, folder)) and folder.endswith('ts')])

"""
# Initialize Dash app
app = Dash(__name__, external_stylesheets=['https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css'])
server = app.server



# Define layout
app.layout = html.Div([
    header(dash_configs['title']),
    scenario_and_tabs(SCENARIOS, plot_files),

    # Add the note at the top about timeslices styled as an alert
    html.Div(
        children=f"Displaying clews model results for {timeslices} timeslices and pypsa model results for hourly resolution",  # Dynamic message
        className="alert alert-info text-center",  # Bootstrap alert class for styling
        style={
            'fontSize': '0.8rem',  # Set the font size
            'color': '#333',  # Set text color for better readability
            'borderRadius': '10px',  # Rounded corners for the alert
            'padding': '2px',  # Padding for better spacing inside the alert
            'boxShadow': '0 2px 10px rgba(0, 0, 0, 0.1)',  # Subtle shadow for a modern touch
        }
    ),

    # Main contents
    html.Div(id='tab-content', className='mt-4')
], className='container')



# Callback to update scenario information
@app.callback(
    Output('scenario-info', 'children'),
    [Input('scenario-dropdown', 'value')]
)
def update_scenario_info(selected_scenario):
    info = info_SCENARIOS.get(selected_scenario, 'No information available for this scenario')
    return html.P(info, className='text-muted')

# Callback to update graph based on tab selection and dropdown value
@app.callback(
    Output('tab-content', 'children'),
    [Input('tabs', 'value'), 
     Input('scenario-dropdown', 'value')] + [Input({'type': 'plot_type-dropdown', 'index': tab_label}, 'value') for tab_label in plot_files]
)
def update_tab_content(selected_tab, selected_scenario, *selected_filenames):
    graphs = []

    if selected_tab == 'resources':
        # Embed the HTML plot as an iframe
        html_file_path = Path("docs/plots/Resource_options.html")
        if html_file_path.exists():
            # Option 1: Display as an iframe
            iframe = html.Iframe(
                srcDoc=html_file_path.read_text(),
                style={
                    'width': '98%',
                    'height': '900px',
                    'border': 'none',
                    'margin': 'auto',
                    'display': 'block',
                    'padding': '10px',
                    'overflow': 'auto',
                }
            )

            # Choose one of the options to return
            return [iframe]
    
    plot_files_list = plot_files[selected_tab]
    selected_filenames = [filename for filename in selected_filenames if filename in plot_files_list]
    if selected_tab == 'demand':
        demand_file_paths = [Path(demand_data) / filename for filename in selected_filenames]
        for file_path, filename in zip(demand_file_paths, selected_filenames):
            if file_path.exists():
                df = pd.read_csv(file_path)
                if all(col in df.columns for col in ["year", "Commercial", "Industrial", "Residential", "Transportation"]):
                    fig=create_demand_plot_simplified(df)
                    graphs.append(dcc.Graph(figure=fig))
                # """ 
                else:
                    aggr_sector_plot,sector_fuels_plot=create_demand_plots(df)

                    # Show figure
                    graphs.append(dcc.Graph(figure=aggr_sector_plot))
                    graphs.append(dcc.Graph(figure=sector_fuels_plot))

                # """   
    if selected_tab == 'dispatch':
            data_2040=pd.read_csv('results/pypsa/generators_t_p_2040.csv',index_col='snapshot',parse_dates=True)
            pypsa_gen_plot=create_generation_plots(data_2040)
            graphs.append(dcc.Graph(figure=pypsa_gen_plot))
    else:
        result_file_paths = [scenario_results_directory / selected_scenario / results_subfolder / filename for filename in selected_filenames]
        for file_path, filename in zip(result_file_paths, selected_filenames):
            if file_path.exists():
                df = pd.read_csv(file_path)
                if filename == "AnnualEmissions.csv":
                    fig = px.line(df, x='YEAR', y=df.columns[-1], title='Emission Trends', markers=False)
                    fig.update_traces(
                        line_color='red',  # Set line color to red
                        opacity=0.9        # Set opacity to 70%
                    )
                    fig.update_xaxes(title_text='Year')
                    fig.update_yaxes(title_text=units_mapping[selected_tab])

                    # Add fixed markers for the BC Emission targets
                    year = [2030, 2040, 2050]
                    emission = [40, 25, 0]  # MTeCO2

                    fig.add_trace(go.Scatter(
                        x=year,
                        y=emission,
                        mode='markers+text',
                        marker=dict(
                            size=12,
                            color='yellow',
                            opacity=0.7,
                            line=dict(width=1, color='orange')
                        ),
                        text=["2030 Target", "2040 Target", "2050 Target"],  # Text labels
                        textposition="top center",
                        name="BC Emission Targets"  # Legend name
                    ))
                    
                    inventory_year=[2021, 2022, 2023]
                    inventory_CO2_data=[63.8,65.6,55.94]# MTeCO2
                    
                    fig.add_trace(go.Scatter(
                        x=inventory_year,
                        y=inventory_CO2_data,
                        mode='markers',
                        marker=dict(
                            size=9,
                            color='blue',
                            opacity=0.5,
                            line=dict(width=1, color='blue')
                        ),
                        # textposition="top center",
                        # text=["2021 actual", "2022 actual", "2023 actual"],  # Text labels
                        name="BC Emission Inventory Report"  # Legend name
                    ))
                    
                    # Ensure full data is shown by setting axis range to fit all data
                    # Update axis titles
                    fig.update_xaxes(
                            title_text='Year',
                            tickmode='linear',  # Show ticks at regular intervals
                            dtick=2,           # Set tick interval (e.g., every 10 years)
                            showgrid=True       # Optional: Show grid lines for better readability
                        )
                    fig.update_yaxes(title_text=units_mapping[selected_tab])
                    graphs.append(dcc.Graph(figure=fig))
                    
                elif filename == "AnnualTechnologyEmission.csv":
                    techs = selected_technologies[selected_tab]
                    years, df_processed = vis_utils.load_and_process_data(file_path, techs)
                    # Create the area plot with negative values preserved
                    fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)

                    # Update x and y axes
                    fig.update_xaxes(title_text='Year')

                    # Set the y-axis range to allow for negative values to be visible
                    fig.update_yaxes(
                        title_text=units_mapping[selected_tab]
                          # Optionally, set to "tozero" if you want the y-axis to start from 0 or leave it dynamic
                    )

                    # Update legend labels (if needed)
                    for tech, label in legend_labels.items():
                        fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())

                    # Append the figure to the list of graphs
                    graphs.append(dcc.Graph(figure=fig))

                else:
                    techs = selected_technologies[selected_tab]
                    years, df_processed = vis_utils.load_and_process_data(file_path, techs)
                    # Ensure years are treated as discrete values (categories)
                    df_processed['Year'] = [str(year) for year in years]  # Convert years to strings for discrete x-axis
                    fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)
                    fig.update_xaxes(title_text='Year')
                    fig.update_yaxes(title_text=units_mapping[selected_tab])
                    for tech, label in legend_labels.items():
                        fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
                    graphs.append(dcc.Graph(figure=fig))
         

    return graphs


# Run the app
if __name__ == '__main__':
    app.run_server(port=8050,debug=True)
