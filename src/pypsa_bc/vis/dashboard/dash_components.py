from dash import html, dcc
from bcnexus import utils

dash_configs: dict = utils.load_config('config/dashboard.yaml')
filenames_mapping = dash_configs['filenames_mapping']

def header(title):
    return html.Div(
        children=[
            html.H1(
                children=title,
                className='text-center my-3 font-weight-bold',
                style={'color': 'darkblue', 'font-size': '1.5rem'}
            ),
            html.Div(
                children=[
                    html.A(
                        "More about BC-Combined-Model",
                        href="https://github.com/DeltaE/BC-CLEWS-Model/wiki",
                        target="_blank",
                        className='btn btn-link btn-sm',
                        style={'position': 'absolute', 'top': '2px', 'right': '2px', 'font-size': '0.9rem', 'color': 'darkblue'}
                    )
                ],
                className='d-flex justify-content-end',
                style={'padding-right': '2px'}
            )
        ],
        style={'background-color': '#f7f7f7', 'padding': '2px 2px', 'position': 'relative', 'margin-bottom': '2px'}
    )


def scenario_and_tabs(SCENARIOS, 
                      plot_file):
    return html.Div(
        children=[
            # Scenario Dropdown
            html.Div(
                children=[
                    html.Label(
                        'Select Scenario',
                        className='font-weight-bold text-primary mb-2',
                        style={'font-size': '0.9rem'}  # Reduced font size
                    ),
                    dcc.Dropdown(
                        id='scenario-dropdown',
                        options=[{'label': scenario, 'value': scenario} for scenario in SCENARIOS],
                        value=SCENARIOS[0],
                        className='form-control',
                        style={
                            'border': 'none',
                            'font-size': '0.5rem',  # Smaller font size
                            'box-shadow': 'none',
                            'padding': '4px 10px',  # Reduced padding
                            'width': '300px'  # Adjusted width
                        },
    
                    ),
                    html.Div(
                        id='scenario-info',
                        className='mt-2',
                        style={'color': '#555', 'font-size': '0.6rem'}
                    )
                ],
                style={'flex': '1', 'margin-right': '10px'}  # Adjust flex behavior and spacing
            ),
            
            # Tabs Component
            html.Div(
                children=[
                    dcc.Tabs(
                        id='tabs',
                        value=list(plot_file.keys())[0],
                        children=[
                            dcc.Tab(
                                label=tab_label,
                                value=tab_label,
                                children=[
                                    html.Div(
                                        children=[
                                            dcc.Dropdown(
                                                id={'type': 'plot_type-dropdown', 'index': tab_label},
                                                options=[{'label': filenames_mapping[filename], 'value': filename} for filename in plot_files],
                                                value=plot_files[0] if plot_files else None,
                                                className='form-control',
                                                style={
                                                    'border': 'none',
                                                    'box-shadow': 'none',
                                                    'font-size': '0.8rem',
                                                    'padding': '2px 4px',
                                                    'height': '10px',
                                                    'line-height': '1.1'
                                                }
                                            )
                                        ],
                                        className='my-1'
                                    )
                                ],
                                style={
                                    'background-color': '#f7f7f7',
                                    'border-radius': '4px',
                                    'padding': '4px 4px',
                                    'margin-bottom': '2px',
                                    'border': 'none'
                                },
                                selected_style={
                                    'background-color': '#007bff',
                                    'color': 'white',
                                    'font-weight': 'bold',
                                    'border-radius': '4px',
                                    'padding': '4px 4px',
                                    'border': 'none'
                                }
                            ) for tab_label, plot_files in plot_file.items()
                        ],
                        className='nav nav-tabs',
                        style={'border': 'none'}
                    )
                ],
                style={'flex': '3'}  # Adjust flex ratio for the tabs
            )
        ],
        style={
            'display': 'flex',  # Flexbox layout
            'align-items': 'flex-start',  # Align items to the top
            'gap': '2px',  # Space between items
            'padding': '2px 2px',
            'margin-bottom': '2px',
            'background-color': '#ffffff'  # Match the background color
        }
    )
