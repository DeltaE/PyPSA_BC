# * Load Packages
import pandas as pd
import plotly.express as px
from pathlib import Path
import plotly.graph_objects as go
import pypsa
import geopandas as gpd
from pypsa_bc import utils

# local packages
# from bc_combined_modelling.vis import vis_utils
from bcnexus.vis.plot_demand import create_demand_plots, create_demand_plot_simplified
from pypsa_bc.vis.plot_pypsa_gen import  create_generation_plots
from pypsa_bc.vis.vis_utils import plot_inter_region_link_usage,visualize_timeseries,load_and_process_data
    
def get_generators_data(pypsa_network_path:str|Path):

    generator_types = ['Discharge', 'RoR', 'Wind', 'Solar', 'Backstop', 'NG', 'biogas', 'biomass','Nuclear']
    generators_timeseries = pd.DataFrame()
    n = pypsa.Network()
    
    if Path(pypsa_network_path).exists:
        pypsa.Network.import_from_netcdf(n=n, path=pypsa_network_path)
    
        for type in generator_types:
            if type in ['Discharge', 'biomass', 'biogas', 'NG']:
                if any(type in col for col in n.links_t.p1.columns):
                    discharge_links = [col for col in n.links_t.p1.columns if type in col]
                    discharge_links_df = n.links_t.p1[discharge_links]
                    discharge_links_aggregates = discharge_links_df.sum(axis=1)
                    generators_timeseries[type] = abs(discharge_links_aggregates)
            else:
                data = n.generators_t.p
                data = data.loc[:, (data != 0).any(axis=0)]
                
                if any(type in col for col in data.columns):
                    columns = [col for col in data.columns if type in col]
                    data_for_type = data[columns]
                    data_type_aggregates = data_for_type.sum(axis=1)
                    generators_timeseries[type] = data_type_aggregates
                    
        generators_timeseries_resampled_D_peak = generators_timeseries.resample('D').max()
        generators_timeseries_resampled_D_avg=generators_timeseries.resample('D').mean()
        return generators_timeseries,generators_timeseries_resampled_D_peak,generators_timeseries_resampled_D_avg
    else:
        utils.print_update(message=f"Network file doesn't exist OR corrupted : {pypsa_network_path}",alert=True)

def get_scenario_plots(Nexus_scenario:str,
                        year:int,
                        timeslices:int,
                        generation_unit:str,
                        consumption_unit:str,
                        create_pypsa_plots:bool) -> dict:
    plots={}
    results_subfolder = f'{timeslices}ts_csvs_gurobi'
    selected_scenario=f'Model_Kotzur_{Nexus_scenario}' # with Storage scenario, to access folders
    nexus_results_root=Path(scenario_results_directory / selected_scenario / results_subfolder)
    vis_save_to_root=Path(f'vis/bccm/{Nexus_scenario}')
    vis_save_to_root.mkdir(parents=True,exist_ok=True)
    
    if create_pypsa_plots:
        pypsa_network_path = f'results/pypsa/BCPypsa_{Nexus_scenario}_{year}.nc'
        
    # PyPSA Plots -----------------------------------------------------------------------------------
        if Path(pypsa_network_path).exists():
            pypsa_plots={}
            plots['pypsa']=pypsa_plots
        # PyPSA Hourly Dispatch locked to 2021
            base_load_year = 2021
            profile = pd.read_csv(f'data/processed_data/load/Hourly_profile_{base_load_year}.csv', index_col='TIME', parse_dates=True)
            fig1 = create_generation_plots(profile,year=base_load_year)
            pypsa_plots["Load_profile"] = fig1 
            fig1.write_html(vis_save_to_root/f'Load_profile_{base_load_year}.html')
            
        # PyPSA Line Usage Plot
            
            line_congestion_plot_save_to=vis_save_to_root/f'Line_Congestion_{year}.png'
            regional_boundaries_GADM_L2=gpd.read_file('data/processed_data/regions/gadm41_Canada_L2_BC.geojson')
            
            n = pypsa.Network()
            pypsa.Network.import_from_netcdf(n=n, path=pypsa_network_path)
        
            data = n.generators_t.p
            data = data.loc[:, (data != 0).any(axis=0)]
            data_sampled = data.resample('D').mean()

            data_filtered_sampled = data_sampled.loc[:, data_sampled.columns.str.contains("Backstop", case=False)]
            data_filtered_sampled.index = data_filtered_sampled.index.map(lambda x: x.replace(year=year))

        # PyPSA Backstops
            pypsa_plots["PyPSA_backstops"] = px.area(data_filtered_sampled, 
                            template="plotly_white", 
                            title=f"Simulated Backstop for Load -{year} | NUC_standard_2035_2xCost_{year}", 
                            labels={"value": "Power", "index": "Time"},
                            color_discrete_sequence=px.colors.sequential.Reds)
            pypsa_plots["PyPSA_backstops"].update_layout(yaxis_title="MW")
            pypsa_plots["PyPSA_backstops"].write_html(vis_save_to_root/f'Backstops_{year}.html')

        # PyPSA Line Congestion Plot
            pypsa_plots["PyPSA_line_congestion"] = plot_inter_region_link_usage(year,
                                                            pypsa_network_path,
                                                            regional_boundaries_GADM_L2,
                                                            line_congestion_plot_save_to,
                                                            plot_title=f'Inter Region Link Usage (Nexus Scenario | {Nexus_scenario})')

        # PyPSA Generation Plots
            generators_timeseries,generators_timeseries_resampled_D_peak,generators_timeseries_resampled_D_avg=get_generators_data(pypsa_network_path)
            pypsa_plots["PyPSA_generation_hourly"]=visualize_timeseries(generators_timeseries,
                                year,
                                'area',
                                f'Simulated Generation {year}| {Nexus_scenario} ',
                                'Load Served (MW)')
            
            pypsa_plots["PyPSA_generation_hourly"].write_html(vis_save_to_root/f'Pypsa_gen_{year}.html')
            pypsa_plots["PyPSA_generation_daily_peak"] =visualize_timeseries(generators_timeseries_resampled_D_peak,
                                year,
                                'area',
                                f'Simulated Generation {year} - Daily Peak| {Nexus_scenario} ',
                                'Load Served (MW)')
            
            pypsa_plots["PyPSA_generation_daily_peak"].write_html(vis_save_to_root/f'Pypsa_gen_{year}_Daily_peaks.html')
            pypsa_plots["PyPSA_generation_daily_average"] =visualize_timeseries(generators_timeseries_resampled_D_avg,
                                year,
                                'area',
                                f'Simulated Generation {year} - Daily Average| {Nexus_scenario} ',
                                'Load Served (MW)')
            
            pypsa_plots["PyPSA_generation_daily_average"].write_html(vis_save_to_root/f'Pypsa_gen_{year}_Daily_mean.html')

        else:
            utils.print_update(message=f"Network file doesn't exist OR corrupted : {Path(pypsa_network_path).name}",alert=True)
            utils.print_update(message=f"Missing/corrupted network file @ : {pypsa_network_path}",alert=True)
    
    else:
        utils.print_update(message=f"Skipping PyPSA plots for scenario :{Nexus_scenario} as per user flag for 'create_pypsa_plots'",alert=True)
    
# Total Demands
    demand_file_names = ['Hourly_profile_2021.csv']
    demand_file_paths = [Path(demand_data) / filename for filename in demand_file_names]
    for file_path, filename in zip(demand_file_paths, demand_file_names):
        if file_path.exists():
            df = pd.read_csv(file_path)
            if all(col in df.columns for col in ["year", "Commercial", "Industrial", "Residential", "Transportation"]):
                
                plots["Sectoral_demand"] = create_demand_plot_simplified(df)
                plots["Sectoral_demand"].write_html(vis_save_to_root/'Sectoral_demand.html')
            else:
                plots["Sectoral_total_demand"], plots["Sectoral_fuel_demand"]  = create_demand_plots(df)
                plots["Sectoral_total_demand"].write_html(vis_save_to_root/'Sectoral_total_demand.html')
                plots["Sectoral_fuel_demand"] .write_html(vis_save_to_root/'Sectoral_fuel_demand.html')

# Nexus plots----------------------------------------------------------------------------------
    if nexus_results_root.exists():
        nexus_plots={}
        nexus_ts_plots={}
        plots['nexus']=nexus_plots
        nexus_plots[f'{timeslices}']=nexus_ts_plots
        
    # Nexus landuse
        land_use_files = [nexus_results_root/ filename for filename in plot_files['landuse']]
        for file_path, filename in zip(land_use_files, plot_files['landuse']):
            df = pd.read_csv(file_path)
            landuse_fig = plot_landuse_for_clusters(df,Nexus_scenario)
            
            nexus_ts_plots["Nexus_landuse"] = landuse_fig 
            landuse_fig.write_html(vis_save_to_root/'Nexus_Landuse.html')
            
    # Nexus Consumption Plots
        consumption_files = [scenario_results_directory / selected_scenario / results_subfolder / filename for filename in plot_files['consumption']]
        plot_unit=units_mapping[f'consumption_{consumption_unit.lower()}']  
        
        for file_path, filename in zip(consumption_files, plot_files['consumption']):
            df = pd.read_csv(file_path)
            if 'gwh' in plot_unit.lower():
                df['VALUE'] = df['VALUE'] * get_PJ_to_GWh_conversion_factor() # Convert PJ to GWh
            fig_sector, fig_fuel = plot_combined_stacked_energy_consumption(df, plot_unit,dash_configs,Nexus_scenario)
            
            nexus_ts_plots["Nexus_sectoral_consumption"] = fig_sector  # Store with a custom name
            fig_sector.write_html(vis_save_to_root/'Nexus_Sectoral_consumption.html')
            
            nexus_ts_plots["Nexus_fuel_consumption"] = fig_fuel  # Store with a custom name
            fig_fuel.write_html(vis_save_to_root/f'Nexus_yearly_Fuel_{plot_unit}.html')

    # Nexus Annual Energy/Generation  Plots
        energy_files = [scenario_results_directory / selected_scenario / results_subfolder / filename for filename in plot_files['energy']]
        techs = selected_technologies['energy']
        
        plot_unit=units_mapping[f'energy_{generation_unit.lower()}']  
        
        for file_path, filename in zip(energy_files, plot_files['energy']):
            plot_prefix=0
            years, df_processed = load_and_process_data(file_path, techs)

            if 'gwh' in plot_unit.lower():
                df_processed= df_processed * get_PJ_to_GWh_conversion_factor() # Convert PJ to GWh
            Nexus_yearly_energy_fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)
            Nexus_yearly_energy_fig.update_xaxes(title_text='Year')
            Nexus_yearly_energy_fig.update_yaxes(title_text=plot_unit)
            Nexus_yearly_energy_fig.update_layout(title_text=f'Energy Generation by Technology ({Nexus_scenario})')
            for tech, label in legend_labels.items():
                Nexus_yearly_energy_fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
                
            nexus_ts_plots[f"Nexus_energy{'' if plot_prefix==0 else plot_prefix}"] = Nexus_yearly_energy_fig
            Nexus_yearly_energy_fig.write_html(vis_save_to_root/f'Nexus_yearly_{plot_unit}.html')
            
            plot_prefix=plot_prefix+1

    # Nexus Capacity Investment  Plots
        result_file_paths = [scenario_results_directory / selected_scenario / results_subfolder / filename for filename in plot_files['capacity']]
        techs = selected_technologies['capacity']
        for file_path, filename in zip(result_file_paths, plot_files['capacity']):
            years, df_processed = load_and_process_data(file_path, techs)
            Nexus_capacity_investmemts_fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)
            Nexus_capacity_investmemts_fig.update_xaxes(title_text='Year')
            Nexus_capacity_investmemts_fig.update_yaxes(title_text=units_mapping['capacity'])
            Nexus_capacity_investmemts_fig.update_layout(title_text=' Total Capacity(' + Nexus_scenario + ')' if 'total' in filename.lower() else 'Capacity Investments (' + Nexus_scenario + ')')
            for tech, label in legend_labels.items():
                Nexus_capacity_investmemts_fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
                
            nexus_ts_plots["Nexus_total_capacity" if 'total' in filename.lower() else "Nexus_capacity_investmemts"] = Nexus_capacity_investmemts_fig
            Nexus_energy_plot_filename="Nexus_total_capacity" if 'total' in filename.lower() else "Nexus_capacity_investmemts"
            Nexus_capacity_investmemts_fig.write_html(vis_save_to_root/f'{Nexus_energy_plot_filename}.html')
            
    # Nexus Generation Timeslices

        RateOfUseByTechnology=pd.read_csv(f'{scenario_results_directory}/Model_Kotzur_{Nexus_scenario}/{results_subfolder}/ProductionByTechnology.csv')
        # Group by 'TECHNOLOGY' and 'TIMESLICE' and calculate the sum of 'VALUE' for each group
        # _grouped_ = RateOfUseByTechnology[RateOfUseByTechnology['TECHNOLOGY'].str.upper().str.startswith('PWR')].groupby(['TECHNOLOGY', 'YEAR', 'TIMESLICE'])['VALUE'].sum().reset_index()
        filtered_RateOfUseByTechnology = RateOfUseByTechnology.loc[
            (RateOfUseByTechnology['TECHNOLOGY'].str.upper().str.startswith('PWR') 
            | 
            RateOfUseByTechnology['TECHNOLOGY'].str.upper().str.startswith('IMPPWR')
            | 
            RateOfUseByTechnology['TECHNOLOGY'].str.upper().str.startswith('EXPPWR'))
            & 
            ~RateOfUseByTechnology['TECHNOLOGY'].str.upper().str.contains('PWRTRN')
        ]
        max_timeslice_digits = len(str(int(timeslices)- 1))
        filtered_RateOfUseByTechnology = filtered_RateOfUseByTechnology.copy()
        filtered_RateOfUseByTechnology.loc[:, 'SEQUENTIAL_TIMESLICE'] = (
            filtered_RateOfUseByTechnology['YEAR'].astype(str) + 
            filtered_RateOfUseByTechnology['TIMESLICE'].astype(str).str.zfill(max_timeslice_digits)
        ).astype(int)
            
        # filtered_RateOfUseByTechnology.loc[:, 'SEQUENTIAL_TS_LABEL'] = filtered_RateOfUseByTechnology['YEAR'].astype(str) + '_TS' + filtered_RateOfUseByTechnology['TIMESLICE'].astype(str)
        filtered_RateOfUseByTechnology.loc[:, 'SEQUENTIAL_TS_LABEL'] = filtered_RateOfUseByTechnology['SEQUENTIAL_TIMESLICE'].astype(str).apply(lambda x: f"{x[:4]} {x[4:]}")
     
        filtered_RateOfUseByTechnology = filtered_RateOfUseByTechnology.sort_values(by='SEQUENTIAL_TIMESLICE')
        filtered_RateOfUseByTechnology['VALUE_MW']=  filtered_RateOfUseByTechnology['VALUE']*get_PJ_to_GWh_conversion_factor(8) #MW= PJ×1E9/seconds in timeslice 

        # Create an area plot for each technology
        Nexus_timeslice_activity_fig = px.area(
            filtered_RateOfUseByTechnology,
            x='SEQUENTIAL_TS_LABEL',
            y='VALUE_MW', #MW
            color='TECHNOLOGY',
            title=f"BCNexus Generation Activity of Power Technologies [{Nexus_scenario}] ({timeslices}ts)",
            labels={'SEQUENTIAL_TS_LABEL': 'Representative Timeslice (Year_TS)', 'VALUE': 'Activity (GW)'},
            color_discrete_map=custom_colors  # Use custom colors if defined
        )

        # Update layout for better visualization
        Nexus_timeslice_activity_fig.update_layout(
            xaxis=dict(tickangle=-90),
            yaxis_title="Activity (MW)",
            legend_title="Technology",
            title_x=0.5,
            template="plotly_white"
        )
        
        nexus_ts_plots["Nexus_timeslice_activity"] = Nexus_timeslice_activity_fig
        # Save the figure as an HTML file
        Nexus_timeslice_activity_fig.write_html(vis_save_to_root/f'Nexus_Generation_activity_{timeslices}_timeslices{Nexus_scenario}.html')

    # Nexus Generation Timeslices Aggregated for Tech
        filtered_RateOfUseByTechnology['TECH_CODE'] = filtered_RateOfUseByTechnology['TECHNOLOGY'].apply(
            lambda x: x[:3] if x.startswith('IMP') or x.startswith('EXP') else (x[3:6] if len(x) > 6 else x)
        )

        RateOfUseByTechnology_grouped_aggr = (
            filtered_RateOfUseByTechnology
            .groupby(['TECH_CODE', 'YEAR', 'SEQUENTIAL_TS_LABEL','SEQUENTIAL_TIMESLICE'], as_index=False)['VALUE_MW']
            .sum()
        )
        RateOfUseByTechnology_grouped_aggr=RateOfUseByTechnology_grouped_aggr.sort_values(by='SEQUENTIAL_TIMESLICE')

        # Create an area plot for each technology
        RateOfUseByTechnology_grouped_aggr_fig = px.area(
            RateOfUseByTechnology_grouped_aggr,
            x='SEQUENTIAL_TS_LABEL',
            y='VALUE_MW', #MW
            color='TECH_CODE',
            title=f"BCNexus Generation Activity of Power Technologies [{Nexus_scenario}] ({timeslices} ts)",
            labels={'SEQUENTIAL_TS_LABEL': 'Representative Timeslice (Year_TS)', 'VALUE': 'Activity (GW)'},
            color_discrete_map=custom_colors  # Use custom colors if defined
        )

        # Map TECH_CODE to more descriptive labels in the legend
        tech_code_labels = {
            'WND': 'Wind',
            'SOL': 'Solar',
            'NUC': 'Nuclear',
            'HYD': 'Hydro',
            'HDG': 'Hydrogen',
            'BIO': 'Biomass',
            'NG': 'Natural Gas',
            'IMP': 'Imported Power',
            'EXP': 'Exported Power',
            # Add more mappings as needed
        }

        for tech_code, label in tech_code_labels.items():
            RateOfUseByTechnology_grouped_aggr_fig.for_each_trace(
            lambda trace: trace.update(name=label) if trace.name == tech_code else ()
            )

        # Update layout for better visualization
        RateOfUseByTechnology_grouped_aggr_fig.update_layout(
            xaxis=dict(tickangle=-90),
            yaxis_title="Activity (MW)",
            legend_title="Technology",
            title_x=0.5,
            template="plotly_white"
        )
        
        nexus_ts_plots["Nexus_timeslice_activity_aggregated"] = RateOfUseByTechnology_grouped_aggr_fig
        # Save the figure as an HTML file
        RateOfUseByTechnology_grouped_aggr_fig.write_html(vis_save_to_root/f'Nexus_Generation_activity_{timeslices}_timeslices_aggregated_tech_{Nexus_scenario}.html')
        
        # Nexus Emission  Plots
        result_file_paths = [scenario_results_directory / selected_scenario / results_subfolder / filename for filename in plot_files['emission']]
        for file_path, filename in zip(result_file_paths, plot_files['emission']):
            if file_path.exists():
                df = pd.read_csv(file_path)
                if filename == "AnnualEmissions.csv":
                    AnnualEmissions_fig = px.line(df, x='YEAR', y=df.columns[-1], title='Emission Trends', markers=False)
                    AnnualEmissions_fig.update_traces(
                        line_color='red',  # Set line color to red
                        opacity=0.9        # Set opacity to 70%
                    )
                    AnnualEmissions_fig.update_xaxes(title_text='Year')
                    AnnualEmissions_fig.update_yaxes(title_text=units_mapping['emission'])

                    # Add fixed markers for the BC Emission targets
                    year = [2030, 2040, 2050]
                    emission = [40, 25, 0]  # MTeCO2

                    AnnualEmissions_fig.add_trace(go.Scatter(
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
                    
                    AnnualEmissions_fig.add_trace(go.Scatter(
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
                    AnnualEmissions_fig.update_xaxes(
                            title_text='Year',
                            tickmode='linear',  # Show ticks at regular intervals
                            dtick=2,           # Set tick interval (e.g., every 10 years)
                            showgrid=True       # Optional: Show grid lines for better readability
                        )
                    AnnualEmissions_fig.update_yaxes(title_text=units_mapping['emission'])
                    
                    nexus_ts_plots["Nexus_emission_trends"] = AnnualEmissions_fig
                    AnnualEmissions_fig.write_html(vis_save_to_root/'Nexus_emission_trends.html')

                elif filename == "AnnualTechnologyEmission.csv":
                    techs = selected_technologies['emission']
                    years, df_processed = load_and_process_data(file_path, techs)
                    # Create the area plot with negative values preserved
                    AnnualEmissions_Techs_fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)

                    # Update x and y axes
                    AnnualEmissions_Techs_fig.update_xaxes(title_text='Year')

                    # Set the y-axis range to allow for negative values to be visible
                    AnnualEmissions_Techs_fig.update_yaxes(
                        title_text=units_mapping['emission']
                            # Optionally, set to "tozero" if you want the y-axis to start from 0 or leave it dynamic
                    )

                    # Update legend labels (if needed)
                    for tech, label in legend_labels.items():
                        AnnualEmissions_Techs_fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
                            
                    nexus_ts_plots["Nexus_emission_technologies"] = AnnualEmissions_Techs_fig
                    AnnualEmissions_Techs_fig.write_html(vis_save_to_root/'Nexus_emission_technologies.html')

                # else:
                #     techs = selected_technologies[emission]
                #     years, df_processed = vis_utils.load_and_process_data(file_path, techs)
                #     fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)
                #     fig.update_xaxes(title_text='Year')
                #     fig.update_yaxes(title_text=units_mapping[emission])
                #     for tech, label in legend_labels.items():
                #         fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
                #         
            
    else:
        utils.print_update(message=f"{nexus_results_root} doesn't exist. Can't generate NEXUS plots",alert=True)

    
    return plots

# tech_key='PWRWND'

# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns


# df=pd.read_csv(f'results/clews/Model_Kotzur_{Nexus_scenario}/8ts_csvs_gurobi/RateOfUseByTechnology.csv')
# df['YEAR'] = df['YEAR'].astype(int)
# df['TIMESLICE'] = df['TIMESLICE'].astype(str)  # for better axis labeling

# # Filter the DataFrame for the desired technology
# filtered_df = df[df['TECHNOLOGY'].str.startswith(tech_key) & df['FUEL'].str.startswith(tech_key)]

# # Check if the filtered DataFrame is not empty
# if not filtered_df.empty:
#     # Create a pivot table
#     df_heatmap = filtered_df.pivot_table(
#         index='YEAR', columns='TIMESLICE', values='VALUE', aggfunc='sum'
#     )

#     # Plot the heatmap
#     plt.figure(figsize=(12, 6))
#     sns.heatmap(df_heatmap, cmap='Greens', annot=True, fmt=".0f")
#     plt.title(f'{tech_key} Generation Heatmap ({Nexus_scenario})')
#     plt.xlabel('Timeslice')
#     # plt.xticks(ticks=range(len(df_heatmap.columns)), labels=df_heatmap.columns, rotation=0)
#     plt.ylabel('Year')
#     plt.show()
# else:
#     print("No data available for the specified technology (PWRWNDB01).")

# techs = dash_configs['technologies']['capacity']
# result_file_paths = [Path(f'results/clews/Model_Kotzur_{Nexus_scenario}/8ts_csvs_gurobi') / filename for filename in plot_files['capacity']]
# result_file_paths=[result_file_paths[0]]
# for file_path, filename in zip(result_file_paths, plot_files['capacity']):
#     years, df_processed = vis_utils.load_and_process_data(file_path, techs)
#     fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)
#     fig.update_xaxes(title_text='Year')
#     fig.update_yaxes(title_text=units_mapping['capacity'])
#     fig.update_layout(title_text=' Total Capacity(' + Nexus_scenario + ')' if 'total' in filename.lower() else 'Capacity Investments (' + Nexus_scenario + ')')
#     for tech, label in legend_labels.items():
#         fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
#     graphs["Nexus_total_capacity" if 'total' in filename.lower() else "Nexus_capacity_investmemts"] = fig
#     Nexus_energy_plot_filename="Nexus_total_capacity" if 'total' in filename.lower() else "Nexus_capacity_investmemts"
#     # fig.write_html(Path('vis')/f'{Nexus_energy_plot_filename}.html')
# fig

# Nexus_scenario=    'NUC_standard_2035_2xCost'
# techs = dash_configs['technologies']['capacity']
# result_file_paths = [Path('results/clews/Model_Kotzur_NUC_standard_2030_2xCost/8ts_csvs_gurobi') / filename for filename in plot_files['capacity']]
# result_file_paths=[result_file_paths[0]]
# for file_path, filename in zip(result_file_paths, plot_files['capacity']):
#     years, df_processed = vis_utils.load_and_process_data(file_path, techs)
#     fig = px.bar(df_processed, x=years, y=techs, title=filenames_mapping[filename], color_discrete_map=custom_colors)
#     fig.update_xaxes(title_text='Year')
#     fig.update_yaxes(title_text=units_mapping['capacity'])
#     fig.update_layout(title_text=' Total Capacity(' + Nexus_scenario + ')' if 'total' in filename.lower() else 'Capacity Investments (' + Nexus_scenario + ')')
#     for tech, label in legend_labels.items():
#         fig.for_each_trace(lambda trace: trace.update(name=label) if trace.name == tech else ())
#     graphs["Nexus_total_capacity" if 'total' in filename.lower() else "Nexus_capacity_investmemts"] = fig
#     Nexus_energy_plot_filename="Nexus_total_capacity" if 'total' in filename.lower() else "Nexus_capacity_investmemts"
#     # fig.write_html(Path('vis')/f'{Nexus_energy_plot_filename}.html')
# fig
