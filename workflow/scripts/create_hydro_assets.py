from bc_power import utils, hydro
import pandas as pd


def get_hydro_data_dict():
    '''
    This function returns the basic data dictionary used as a template for storing
    the features of each hydro generation asset.
    
    '''
    # data_dict: holds the features which will be required in the final hydro_assets.csv
    data_dict = { # replace later with reading a configuration file with the targetted columns
                "asset_id":"DEFAULT",
                "connecting_node_code":"DEFAULT",
                "num_of_units":"DEFAULT",
                "latitude":"DEFAULT",
                "longitude":"DEFAULT",
                "capacity":"DEFAULT",
                "capacity_factor":"DEFAULT",
                "annual_avg_energy":"DEFAULT",
                "p_min":"DEFAULT",
                "nominal_head":"DEFAULT",
                "max_spill":"DEFAULT",
                "max_water_discharge":"DEFAULT",
                "min_water_discharge":"DEFAULT",
                "gen_type":"DEFAULT",
                "cascade_group":"DEFAULT",
                "cascade_order":"DEFAULT",
                "upper_reservoir_id":"DEFAULT",
                "lower_reservoir_id":"DEFAULT",
                "max_level":"DEFAULT",
                "min_level":"DEFAULT",
                "min_storage":"DEFAULT",
                "max_storage":"DEFAULT",
                "ramp_up":"DEFAULT",
                "ramp_down":"DEFAULT",
                "min_up":"DEFAULT",
                "min_down":"DEFAULT",
                "min_gen":"DEFAULT",
                "spinning_reserve_capability":"DEFAULT",
                "forced_outage_rate":"DEFAULT",
                "planned_outage_rate":"DEFAULT",
                "start_up_cost_cold":"DEFAULT",
                "shutdown_cost":"DEFAULT",
                "typical_plant_size_MW":"DEFAULT",
                "capital_cost_USD_per_kW":"DEFAULT",
                "service_life_years":"DEFAULT",
                "fixed_om_cost_USD_per_MWyear":"DEFAULT",
                "variable_om_cost_USD_per_MWh":"DEFAULT",
                "average_fuel_price_USD_per_MMBtu":"DEFAULT",
                "carbon_emissions_tCO2eq_per_MWh":"DEFAULT",
                "start_year":"DEFAULT",
                "closure_year":"DEFAULT"}

    return data_dict


def add_generators_features(row, cid_dict):
    '''
    This function extracts generator features from the generators.csv table
    from CODERS to dictionary which is merging hydro relevant data to create a
    hydro assets file.
    '''
    # asset_id conversion

    cid_dict["asset_id"] = utils.convert_cid_2_aid(row['gen_node_code'], row["connecting_node_code"])
    cid_dict["connecting_node_code"] = row["connecting_node_code"]
    cid_dict["num_of_units"] = row["total_num_of_units"]
    cid_dict["latitude"] = row["latitude"]
    cid_dict["longitude"] = row["longitude"]
    cid_dict["capacity"] = row["install_capacity_in_mw"] # Could potentially used effective instead of installed capacity
    cid_dict["annual_avg_energy"] = row["annual_avg_energy_unit_in_gwh/y"]
    cid_dict["gen_type"] = row["gen_type"].lower()
    cid_dict["start_year"] = row["start_year"]
    cid_dict["closure_year"] = row["closure_year"]


def add_cascade_features(row, cid_dict):
    '''
    This function extracts generator features from the hydro_cascade.csv table
    from CODERS to dictionary which is merging hydro relevant data to create a
    hydro assets file.
    '''
    cid_dict["nominal_head"] = "DEFAULT" # Too sparse
    cid_dict["max_spill"] = float(row["max_spill"]) # m3/s
    cid_dict["max_water_discharge"] = row["max_water_discharge"] # m3/s
    cid_dict["min_water_discharge"] = row["min_water_discharge"] # m3/s
    cid_dict["cascade_group"] = row["cascade_group_name"]
    cid_dict["cascade_order"] = row["number"]
    cid_dict["max_storage"] = row["max_storage"]
    cid_dict["min_storage"] = row["min_storage"]

def add_hydro_existing_features(row, cid_dict):
    '''
    This function extracts generator features from the existing.csv table
    from CODERS to dictionary which is merging hydro relevant data to create a
    hydro assets file.
    '''
    cid_dict["upper_reservoir_id"] = row["upper_storage_name"]
    cid_dict["lower_reservoir_id"] = row["lower_storage_name"]
    cid_dict["max_level"] = row["maximum_level"]
    cid_dict["min_level"] = row["minimum_level"]

def add_gen_generic_features(row, cid_dict):
    '''
    This function extracts generator features from the generation_generic.csv table
    from CODERS to dictionary which is merging hydro relevant data to create a
    hydro assets file.
    UPDATE: Likely, may use try other sources from ramping constraints and etc.. in future. This function will need to be updated
            accordingly to accommodate those changes.
    '''
    cid_dict["ramp_up"] = cid_dict["capacity"] * row["ramp_rate_percent_per_min"] # MW/min
    cid_dict["ramp_down"] = cid_dict["capacity"] * row["ramp_rate_percent_per_min"] # MW/min (symmetric assume..)
    cid_dict["min_up"] = row["min_up_time_hours"] # hrs
    cid_dict["min_down"] = row["min_down_time_hours"] # hrs
    cid_dict["min_gen"] = row["min_plant_load"] * cid_dict["capacity"] # hrs
    cid_dict["spinning_reserve_capability"] = row["spinning_reserve_capability"] # ??? (Needs to be corrected...)
    cid_dict["forced_outage_rate"] = row["forced_outage_rate"]
    cid_dict["planned_outage_rate"] = row["planned_outage_rate"]
    cid_dict["start_up_cost_cold"] = row["cold_start_up_costs_$_per_Mwcap"]
    cid_dict["shutdown_cost"] = row["shutdown_cost"]
    cid_dict["typical_plant_size_MW"] = row["typical_plant_size_MW"]
    cid_dict["capital_cost_USD_per_kW"] = row["capital_cost_USD_per_kW"]
    cid_dict["service_life_years"] = row["service_life_years"]
    cid_dict["fixed_om_cost_USD_per_MWyear"] = row["fixed_om_cost_USD_per_MWyear"]
    cid_dict["variable_om_cost_USD_per_MWh"] = row["variable_om_cost_USD_per_MWh"]
    cid_dict["average_fuel_price_USD_per_MMBtu"] = row["average_fuel_price_USD_per_MMBtu"]
    cid_dict["carbon_emissions_tCO2eq_per_MWh"] = row["carbon_emissions_tCO2eq_per_MWh"]

def get_features_generators(generators, component_dict, data_dict, hydro_types):
    '''
    This function gets all features from the CODERS generators.csv file for only hydro assets.
    '''
    mask = (generators['province'] == "BC") & (generators['gen_type'].apply(lambda x: x.lower() in hydro_types))
    generators[mask]
    for idx,row in generators[mask].iterrows():
        cid = row["gen_node_code"]
        if cid in component_dict:
            print("ERROR Non-unique: Component ID already added.")
            break
        component_dict[cid] = data_dict.copy()
        add_generators_features(row, component_dict[cid])

def get_feature_cascade(cascade_data, component_dict):
    '''
    This function gets all features from the CODERS hydro_cascade.csv
    '''
    for idx,row in cascade_data.iterrows():
        cid = row["gen_node_code"]
        if cid not in component_dict:
            print("ERROR Non-unique: Component ID already added.")
            break
        add_cascade_features(row, component_dict[cid])

def get_feature_existing(hydro_e_data, component_dict):
    '''
    This function gets all features from the CODERS hydro_existing.csv
    '''
    mask = (hydro_e_data['Province'] == "BC") 
    for idx,row in hydro_e_data[mask].iterrows():
        for cid in component_dict.keys(): # many-to-one
            if cid not in component_dict:
                print("ERROR Non-unique: Component ID already added.")
                break
            if component_dict[cid]["asset_id"] == row["connecting_node_code"]:
                add_hydro_existing_features(row, component_dict[cid])

def get_feature_generic(gen_generic,component_dict):
    '''
    This function gets all features from the CODERS generation_generic.csv
    '''
    gen_generic.dropna(how="all",inplace=True)
    for cid,features in component_dict.items():
        for idx,row in gen_generic.iterrows():
            if features["gen_type"] == row["generation_type"].lower():
                add_gen_generic_features(row, features)

def create_df_hydro_gen(component_dict):
    '''
    This function creates a dataframe from the component dictionary. This dictionary contains
    all the information of the hydro generation assets from CODERs.
    Note: This function fixes a few edge case asset ids.
    Note: This function show eventually have the custom entries pulled from hard coded in another location.
    '''
    # unwrap and create csv 
    # and Custom fix for CID: BC_USR00_GEN, BC_USR02_GEN
    df_hydro_gen = pd.DataFrame.from_dict(component_dict,orient='index')
    df_hydro_gen.reset_index(inplace=True)
    df_hydro_gen.rename({"index":"component_id"}, inplace=True, axis=1)
    cid_2_aid_dict = {"BC_USR00_GEN":"BC_USR_GSS","BC_USR02_GEN":"BC_LMN_GSS"} # custom fix
    for cid,aid in cid_2_aid_dict.items():
        ind = df_hydro_gen[df_hydro_gen["component_id"] == cid].index[0]
        df_hydro_gen.at[ind,"asset_id"] = aid

    return df_hydro_gen


def add_hydro_type(hydro_sites, res_wup_data, cfg):
    '''
    This function is used to add a column type to the hydro generation dataframe.
    This column can then be used to split the hydro generatio assets into ones needing
    inflow series calculations and those needing power availability series calculated.
    Currently the following situations are supported for each type:
    hydro_sites: Dataframe of hydroelectric generation assets.
    res_wup_data: Dataframe of hydroelectric reservoir assets.
    Inflow Series (Reservoir type):
    i) Reservoirs with WUP statistics.
    ii) Reservoirs downstream in a cascade w/o WUP statistics.
    Availability Series (RoR type):
    i) RoR assets.
    ii) Reservoirs w/o WUP statistics and are not downstream of other modelled reservoirs.
    Tags of "reservoir" and "ror" will be used respectively
    '''
    
    new_hydro_col = "hydro_type"

    # a) Split into cascade and non-cascade
    rid_list = res_wup_data["asset_id"].tolist()
    mask_cascade = hydro_sites['upper_reservoir_id'].apply(lambda x: x in rid_list)
    cascade_sites = hydro_sites[mask_cascade]
    ror_index = hydro_sites[~mask_cascade].index
    hydro_sites.loc[ror_index,new_hydro_col] = "ror"

    # b) Get mask of cascaded reservoirs with WUP statistics
    mask_res_wup = cascade_sites['upper_reservoir_id'].apply(lambda x: hydro.check_wup_exists(x,
                                        cfg['bc_hydro']['inflow_tables']))
    mask_res_wup_index = cascade_sites[mask_res_wup].index

    hydro_sites.loc[mask_res_wup_index,new_hydro_col] = "reservoir"

    # c) Items in cascade_sites[~mask_res_wup] need to be checked for upstream or downstream
    # Upstream head reservoirs (no reservoirs upstream of these) with no WUP statistics
    mask_res_up_no_wup = cascade_sites['upper_reservoir_id'].apply(lambda x: 
                        hydro.check_head_reservoir(x, cascade_sites)) & (~mask_res_wup)
    mask_res_up_no_wup_index = cascade_sites[mask_res_up_no_wup].index # BEWARE

    hydro_sites.loc[mask_res_up_no_wup_index,new_hydro_col] = "ror"

    # Downstream reservoirs with no WUP statistics
    mask_res_down_no_wup = (~cascade_sites['upper_reservoir_id'].apply(lambda x: 
                        hydro.check_head_reservoir(x, cascade_sites))) & (~mask_res_wup)
    mask_res_down_no_wup_index = cascade_sites[mask_res_down_no_wup].index

    hydro_sites.loc[mask_res_down_no_wup_index,new_hydro_col] = "reservoir-impute"

    return hydro_sites

def main():
    '''
    Description: This script is for creating a single csv file with all
    the needed technical information regarding the hydro assets.     
    
    '''
    # Description: This script is for creating a single csv file with all
    #  the needed technical information regarding the hydro assets. 


    # Read in configuration file
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    # write path + file
    df_hydro_path = cfg['hydro_prep']["hydro_generation"]
    df_res_path = cfg['hydro_prep']["hydro_reservoir"]

    # Read in files
    gen_generic = pd.read_csv(cfg["coders"]["gen_generic"])
    generators = pd.read_csv(cfg["coders"]["generators"])
    hydro_e_data = pd.read_csv(cfg["coders"]["hydro_existing"])
    cascade_data = pd.read_csv(cfg["coders"]["hydro_cascade"])
    gen_wup_data= pd.read_csv(cfg["bc_hydro"]["generation_wup"])
    res_wup_data= pd.read_csv(cfg["bc_hydro"]["reservoir_wup"])

    # get templates
    data_dict = get_hydro_data_dict()
    hydro_types = set(["hydro_daily","hydro_run","hydro_monthly","hydro_annual"])
    component_dict = {}

    # (i) get features from generators.csv
    get_features_generators(generators, component_dict, data_dict, hydro_types)

    # (ii) get features from hydro_cascade.csv
    get_feature_cascade(cascade_data, component_dict)

    # (iii) get features from hydro existing
    get_feature_existing(hydro_e_data, component_dict)

    # (iv) generic hydro dataset
    get_feature_generic(gen_generic,component_dict)

    # (v)Create df and customzied imputation of select areas
    df_hydro_gen = create_df_hydro_gen(component_dict)

    # (vi) 
    # a) update hydro technical parameters based on Water Use Plan (WUP) generation data
    # b) create csv of hydro generation assets
    for ind,row in gen_wup_data.iterrows():
        temp_mask = df_hydro_gen["asset_id"] == row["asset_id"]
        for ser_ind,ser_val in row[row.notnull()].items():
            df_hydro_gen.loc[temp_mask,ser_ind] = ser_val # Use non-empty values from the generation WUP extacted data
    
    # (vii)
    # a) code determines what type of time-series and modelling structure is required for each asset
    # Type 1) Hydro reservoirs apart of cascades with WUP data or those w/o WUP data and are downstream. (Need Inflow)
    # Type 2) RoR or smaller reservoirs which feed into cascades (these do not have WUP data) (Need power availability)
    add_hydro_type(df_hydro_gen, res_wup_data, cfg)

    # Write files
    df_hydro_gen.to_csv(df_hydro_path, index=False)
    res_wup_data.to_csv(df_res_path, index=False)

    
 
if __name__ == '__main__':
    main()