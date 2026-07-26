from pypsa_bc import utils,hydro
import pandas as pd
from pathlib import Path
import warnings
import sys

# handles the config loading centrally
from pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()

# Suppress all warnings
warnings.filterwarnings("ignore")

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
                # "typical_plant_size_MW":"DEFAULT", # deprecated in CODERS, EL_20270724
                "capital_cost_CAD_per_kW":"DEFAULT",
                "service_life_years":"DEFAULT",
                "fixed_om_costs":"DEFAULT",
                "variable_om_costs":"DEFAULT",
                "average_fuel_price_CAD_per_MMBtu":"DEFAULT",
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

    cid_dict["asset_id"] = utils.convert_cid_2_aid(row['generation_unit_code'], row["network_node_code"])
    cid_dict["connecting_node_code"] = row["network_node_code"]
    cid_dict["num_of_units"] = row["total_facility_generation_units"]
    cid_dict["latitude"] = row["latitude"]
    cid_dict["longitude"] = row["longitude"]
    cid_dict["capacity"] = row["unit_installed_capacity"] # facility_installed_capacity instead of unit_installed_capacity
    cid_dict["annual_avg_energy"] = row["facility_average_annual_energy"]
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
    # Defensive: the current CODERS hydro_existing may lack these columns
    # (reservoir topology + levels now come from the WUP files). Set only what exists.
    if "upper_storage_name" in row:
        cid_dict["upper_reservoir_id"] = row["upper_storage_name"]
    if "lower_storage_name" in row:
        cid_dict["lower_reservoir_id"] = row["lower_storage_name"]
    if "maximum_level" in row:
        cid_dict["max_level"] = row["maximum_level"]
    if "minimum_level" in row:
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
    cid_dict["start_up_cost_cold"] = row["startup_cost"]
    cid_dict["shutdown_cost"] = row["shutdown_cost"]
    # cid_dict["typical_plant_size_MW"] = row["typical_plant_size_MW"] # deprecated in CODERS, EL_20270724
    cid_dict["capital_cost_CAD_per_kW"] = row["capital_cost_CAD_per_kW"]
    cid_dict["service_life_years"] = row["service_life"]
    cid_dict["fixed_om_costs"] = row["fixed_om_costs"]
    cid_dict["variable_om_costs"] = row["variable_om_costs"]
    cid_dict["average_fuel_price_CAD_per_MMBtu"] = row["average_fuel_price_CAD_per_MMBtu"]
    cid_dict["carbon_emissions_tCO2eq_per_MWh"] = row["carbon_emissions"]

def get_features_generators(generators, component_dict, data_dict, hydro_types, cfg=None):
    '''
    This function gets all features from the CODERS generators.csv file for only hydro assets.
    '''
    # mask = (generators['province'].apply(lambda x: x in cfg['output']['prepare_base_network']['regions'])) & (generators['gen_type'].apply(lambda x: x.lower() in hydro_types))
    mask = generators['gen_type'].apply(lambda x: x.lower() in hydro_types)
    
    for idx,row in generators[mask].iterrows():
        cid = row["generation_unit_code"]
        if cid in component_dict:
            print("ERROR Non-unique: Component ID {} already added".format(cid))
            break
        component_dict[cid] = data_dict.copy()
        add_generators_features(row, component_dict[cid])

def get_feature_cascade(cascade_data, component_dict):
    '''
    This function gets all features from the CODERS hydro_cascade.csv
    '''
    for idx,row in cascade_data.iterrows():
        cid = row["gen_node_code"] # Modified: 2024-10-01 for CODERS naming convention update. Hydro tables were not updated and have old naming still.
        if cid not in component_dict:
            print(cid)
            print("ERROR Non-unique: Component ID already added.")
            break
        add_cascade_features(row, component_dict[cid])

def get_feature_existing(hydro_e_data, component_dict, cfg=None):
    '''
    This function gets all features from the CODERS hydro_existing.csv
    '''
    # mask = (hydro_e_data['Province'].apply(lambda x: x in cfg['output']['prepare_base_network']['regions'])) 
    for idx,row in hydro_e_data.iterrows():
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
            if features["gen_type"] == row["gen_type"].lower():
                add_gen_generic_features(row, features)

def create_df_hydro_gen(component_dict, cfg):
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
    
    if "BC" in pypsa_aparser.base_network_cfg['regions']:  # regions moved to base_network.yaml
        cid_2_aid_dict = {"BC_USR00_GEN":"BC_USR_GSS","BC_USR02_GEN":"BC_LMN_GSS"} # custom fix
        for cid,aid in cid_2_aid_dict.items():
            ind = df_hydro_gen[df_hydro_gen["component_id"] == cid].index[0]
            df_hydro_gen.at[ind,"asset_id"] = aid

    return df_hydro_gen


def add_hydro_type(hydro_sites, res_wup_data, inflow_tables):
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
                                        inflow_tables))
    mask_res_wup_index = cascade_sites[mask_res_wup].index

    hydro_sites.loc[mask_res_wup_index,new_hydro_col] = "reservoir"

    # c) Items in cascade_sites[~mask_res_wup] need to be checked for upstream or downstream
    # i) Upstream head reservoirs (no reservoirs upstream of these) with no WUP statistics
    mask_res_up_no_wup = cascade_sites['upper_reservoir_id'].apply(lambda x: 
                        hydro.check_head_reservoir(x, cascade_sites)) & (~mask_res_wup)
    mask_res_up_no_wup_index = cascade_sites[mask_res_up_no_wup].index # BEWARE

    hydro_sites.loc[mask_res_up_no_wup_index,new_hydro_col] = "ror"

    # ii) Downstream reservoirs with no WUP statistics
    mask_res_down_no_wup = (~cascade_sites['upper_reservoir_id'].apply(lambda x: 
                        hydro.check_head_reservoir(x, cascade_sites))) & (~mask_res_wup)
    mask_res_down_no_wup_index = cascade_sites[mask_res_down_no_wup].index

    hydro_sites.loc[mask_res_down_no_wup_index,new_hydro_col] = "reservoir-impute"

    # d) Custom sites which require potential ror inflow series but also
    # will be modelled in such a way to account for their water discharges.
    # WHN and WDN both need modifications
    asset_list = ['WHN','WDN']
    for aid in asset_list:
        aid_mask = hydro_sites["asset_id"].apply(lambda x: x.split("_")[1] == aid)
        aid_mask_index = hydro_sites[aid_mask].index
        hydro_sites.loc[aid_mask_index, new_hydro_col] = "ror-water"

    # Waneta (WAN) will be modelled as a RoR asset (Same as the expansion which is already captured in above code.)
    aid_mask = hydro_sites["asset_id"].apply(lambda x: x.split("_")[1] == "WAN")
    aid_mask_index = hydro_sites[aid_mask].index
    hydro_sites.loc[aid_mask_index, new_hydro_col] = "ror"

    return hydro_sites

def custom_bridge_agg(df):
    '''
    Aggregate 2 parrellel generation stations "BC_BR0_DFS" and "BC_BR0_GSS",
    these share the same upper and lower reservoir, therefore, from the modelling perspective our parrallel turbines.
    
    NOTE: The hydro_cascade.csv file and the generators.csv file have not both been updated simulataneously.
          This causes an issue with the connection node codes not matching for the same assets between the 2 files.
    '''
    # 1) Identify (skip gracefully if the bridge pair isn't present)
    if 'BC_BR1_GSS' not in df.index or 'BC_BR2_GSS' not in df.index:
        return df
    update_idx = df[df.index == 'BC_BR1_GSS'].index[0]
    remove_idx = df[df.index == 'BC_BR2_GSS'].index[0]

    # 2) Aggregate
    df.loc[update_idx,'capacity'] += df.loc[remove_idx,'capacity']
    df.loc[update_idx,'annual_avg_energy'] += df.loc[remove_idx,'annual_avg_energy']
    df.loc[update_idx,'max_water_discharge'] = float(df.loc[update_idx,'max_water_discharge']) + float(df.loc[remove_idx,'max_water_discharge'])
    df.loc[update_idx,'max_spill'] =  float(df.loc[update_idx,'max_spill']) + float(df.loc[remove_idx,'max_spill'])
    df.loc[update_idx,'num_of_units'] += df.loc[remove_idx,'num_of_units']
    df.loc[update_idx,'ramp_up'] += df.loc[remove_idx,'ramp_up']
    df.loc[update_idx,'ramp_down'] += df.loc[remove_idx,'ramp_down']

    # 3) Delete old entries
    df.drop(index=remove_idx, inplace=True)

    return df


def main():
    '''
    Description: This script is for creating a single csv file with all
    the needed technical information regarding the hydro assets.     
    
    '''
    # Description: This script is for creating a single csv file with all
    #  the needed technical information regarding the hydro assets. 

    utils.print_update(level=1,message="Preparing hydro assets and reservoir data for PyPSA-BC...")
    
    # Read in configuration file
    # config_file = r"config/data.yaml"
    cfg = pypsa_aparser.data_cfg 
    
    
    # write path + file
    df_hydro_path = cfg["output"]["create_hydro_assets"]["hydro_generation"]
    df_res_path = cfg["output"]["create_hydro_assets"]["hydro_reservoir"]

    # Read CODERS via the shared client (coders.yaml aliasing applies)
    from pypsa_bc.data.coders import get_coders
    coders = get_coders()
    gen_generic = coders.load_table("generation_generic", as_gdf=False)
    generators = coders.load_table("generators", as_gdf=False)
    # hydro_existing = hydro subset of generators (from the coders module, aliased);
    # hydro_cascade  = reconstructed table (produced by prepare_hydro via to_legacy_cascade).
    hydro_e_data = coders.existing_hydro
    cascade_data = pd.read_csv(cfg["output"]["create_hydro_assets"]["hydro_cascade"])
    
    # Modified: 2024-10-01, Fix hydro cascade files to ensure naming matches coders updated asset names...
    # Likely maintainers of CODERS will fix this issue in the future, however, manually fixed for time being.
    # map of old to new gen_node_codes for the hydro cascade file
    # Modification A: Cascade hydro file
    bridge_codes = {'BC_BR0101_GEN':'BC_BR101_GEN','BC_BR0102_GEN':'BC_BR102_GEN','BC_BR0103_GEN':'BC_BR103_GEN','BC_BR0104_GEN':'BC_BR104_GEN',
                   'BC_BR0201_GEN':'BC_BR201_GEN','BC_BR0202_GEN':'BC_BR202_GEN','BC_BR0203_GEN':'BC_BR203_GEN','BC_BR0204_GEN':'BC_BR204_GEN'}
    mask = cascade_data['gen_node_code'].isin(bridge_codes)
    cascade_data.loc[mask,'gen_node_code'] = cascade_data.loc[mask,'gen_node_code'].apply(lambda old_cold: bridge_codes[old_cold])

    # Modification B: Existing hydro file
    bridge_codes = {'BC_BR1_DFS':'BC_BR1_GSS'}
    mask = cascade_data['connecting_node_code'].isin(bridge_codes)
    cascade_data.loc[mask,'connecting_node_code'] = cascade_data.loc[mask,'connecting_node_code'].apply(lambda old_cold: bridge_codes[old_cold])
    

    # get templates
    data_dict = get_hydro_data_dict()
    hydro_types = set(["hydro_daily","hydro_run","hydro_monthly","hydro_annual"])
    component_dict = {}

    # (i) get features from generators.csv
    get_features_generators(generators, component_dict, data_dict, hydro_types, cfg)

    # (iii) get features from hydro existing
    get_feature_existing(hydro_e_data, component_dict, cfg)

    # (iv) generic hydro dataset
    get_feature_generic(gen_generic,component_dict)

    # (v) Create df and customzied imputation of select areas
    df_hydro_gen = create_df_hydro_gen(component_dict, cfg)
    df_hydro_gen["hydro_type"] = 'ror' # NOTE: Hydro types all 'ror' by default

    # (ii) get features from hydro_cascade.csv
    get_feature_cascade(cascade_data, component_dict)
    
    gen_wup_data = pd.read_csv(cfg["inventory"]["gen_wup"]) # EL: source ???
    res_wup_data = pd.read_csv(cfg["inventory"]["res_wup"]) # EL: source ???
    inflow_tables = cfg["inventory"]["inflow_tables"] # EL: source ???
    # (vi) 
    # a) update hydro technical parameters based on Water Use Plan (WUP) generation data
    # b) create csv of hydro generation assets
    for ind,row in gen_wup_data.iterrows():
        temp_mask = df_hydro_gen["asset_id"] == row["asset_id"]
        for ser_ind,ser_val in row[row.notnull()].items():
            if ser_ind not in df_hydro_gen.columns:
                continue  # skip provenance/junk WUP columns (sources, Notes, Unnamed:*)
            # Template columns start as strings ("DEFAULT"); allow numeric WUP overrides
            # (pandas 2.x StringDtype rejects setting a float into a string column).
            if df_hydro_gen[ser_ind].dtype != object:
                df_hydro_gen[ser_ind] = df_hydro_gen[ser_ind].astype(object)
            df_hydro_gen.loc[temp_mask, ser_ind] = ser_val  # non-empty WUP overrides
    
    # (vii)
    # a) code determines what type of time-series and modelling structure is required for each asset
    # Type 1) Hydro reservoirs apart of cascades with WUP data or those w/o WUP data and are downstream. (Need Inflow)
    # Type 2) RoR or smaller reservoirs which feed into cascades (these do not have WUP data) (Need power availability)
    add_hydro_type(df_hydro_gen, res_wup_data, inflow_tables)

    # BC only implemented so far for this....
    from pathlib import Path as _Path
    _Path(df_res_path).parent.mkdir(parents=True, exist_ok=True)
    res_wup_data.to_csv(df_res_path, index=False)
    utils.print_update(level=2,message=f"reservoir data saved to: {df_res_path}")


    # (viii)
    # Aggregate hydro turbines into singular assets (sum capacity, first for the rest).
    # groupby.agg is robust to the pandas 2.x change where groupby.apply drops the key.
    sum_cols = {"capacity": "sum"}
    other_cols = {c: "first" for c in df_hydro_gen.columns
                  if c not in sum_cols and c != "asset_id"}
    df_hydro_gen_final = (df_hydro_gen.groupby("asset_id", as_index=False)
                          .agg({**sum_cols, **other_cols})).set_index("asset_id")

    # (ix)
    # Custom aggregation of parallel generation assets in the bridge cascade
    df_hydro_gen_final = custom_bridge_agg(df_hydro_gen_final)

    # Write files
    from pathlib import Path as _Path
    _Path(df_hydro_path).parent.mkdir(parents=True, exist_ok=True)
    df_hydro_gen_final.to_csv(df_hydro_path)
    utils.print_update(level=2,message=f"hydro generator data saved to: {df_hydro_path}")
 
if __name__ == '__main__':
    main()