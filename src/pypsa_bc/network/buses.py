import math
import pandas as pd
from pypsa_bc import utils
from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.reporting.assumptions import log_assumption
from pypsa_bc.reporting.logger import get_logger
from pathlib import Path
from pypsa_bc.network import lines

log = get_logger("buses")

# instance of Attributes Parser to access the base_network config and other attributes
aparser=AttributesParser()
BASE_NETWORK_CFG=aparser.base_network_cfg

# shared CODERSData instance (one per run; tables cached in-memory)
coders_data=get_coders()
REQUIRED_COLUMNS = coders_data.required.get('substations', [])
META = {"module":  Path(__file__).stem , 
        "scenario": aparser.get_scenario,
        "path": aparser.assumption_report_save_to}

INTERTIE_SUFFIXES = {"IPT", "INT"}       # AB ties (IPT) and US/BPA ties (INT)

def node_role(code: str) -> str:
    return "intertie" if code.split("_")[-1] in INTERTIE_SUFFIXES else "internal"



def check_missing_buses(prepared_substations: pd.DataFrame, 
                        prepared_lines: pd.DataFrame):
    '''
    Checks for substations which are missing from the lines dataset.
    The only bus like this for BC is BC_WAX_GSS which should be connected to BC_WAN.
    df_buses_bc: DF being prepared for saving and loading into PyPSA.
    df_sub_bc: DF from CODERS of all substations.
    '''
    # Checking to BC_WAX_GSS is the only one missing
    sub_unique_codes = prepared_substations["node_code"].apply(lambda x: x.split('_')[1].strip(' ')).unique().tolist()
    line_unique_codes_start = prepared_lines["starting_node_code"].apply(lambda x: x.split('_')[1].strip(' ')).unique().tolist()
    line_unique_codes_end = prepared_lines["ending_node_code"].apply(lambda x: x.split('_')[1].strip(' ')).unique().tolist()
    line_unique_codes = list(set(line_unique_codes_start + line_unique_codes_end))

    for code in sub_unique_codes:
        if code not in line_unique_codes:
            log.warning(f"substation code {code} is missing from the lines dataset")


def create_bus_df(prepared_lines: pd.DataFrame, 
                  prepared_substations: pd.DataFrame, 
                  prepared_generators: pd.DataFrame):
    '''
    This function will create an initial DataFrame of buses for PyPSA_BC from a DataFrame of lines.
    When creating the buses it
    The lines DF contains the node names for the buses, nominal voltage, and the carrier is implicitly added.
    '''
    # name = []
    # x = []
    # y = []
    # type = []
    # v_nom = []
    data_dict = {}

    # (1) Add buses based on line nodes
    for _,line in prepared_lines.iterrows():
        # Search for match between line and substation
        for node_code in [line["starting_node_code"], line["ending_node_code"]]:
            # --- guard: intertie boundary nodes have no gen/sub and no coords ---
            if node_code.split("_")[-1] in INTERTIE_SUFFIXES:
                # tag = node_code.split("_")[1][:4]              # 'ABBC03' -> 'ABBC'
                # bus = INTERTIE_BUS.get(tag, "AB")              # AB / WA boundary bus
                # x, y = INTERTIE_COORD[bus]                     # fixed placeholder coord
                
                log_assumption(parameter="AB-BC intertie handling", 
                           value="terminated on AB boundary bus; no import path",
                           rationale="BC-internal scope; interties out of scope for this run",
                           source="lines end node BC_ABBC03_IPT",
                           **META)
                pass
                # return bus, x, y
            else:
                bus_name, bus_x, bus_y = get_bus_name_x_y(line, node_code, prepared_substations, prepared_generators)
                if bus_name not in data_dict: # Avoid duplication (Change to dictionary)
                    data_dict[bus_name] = {'x':bus_x, 'y':bus_y, 'type':line['type'], 'v_nom':line['v_nom']}
                    # name.append(bus_name) # i.e. 230_AAL
                    # x.append(bus_x)
                    # y.append(bus_y)
                    # type.append(line['type'])
                    # v_nom.append(line['v_nom'])
    
    df_buses = pd.DataFrame.from_dict(data_dict,orient='index').reset_index().rename(columns={'index':'name'})
    # df_buses = pd.DataFrame()
    # df_buses = pd.DataFrame()
    # df_buses['name'] = name
    # df_buses['x'] = x
    # df_buses['y'] = y
    # df_buses['type'] = type
    # df_buses['v_nom'] = v_nom

    return df_buses

def get_bus_name_x_y(line, node_code, df_substations, generators):
    '''
    This function finds the correct name for a bus from the node dataset
    line: Line row from CODERS lines dataframe.
    node_code: Name of start/end node from CODERS lines dataframe.
    df_substations: Substations dataframe from CODERS.
    return bus_name: Unique name to use for the bus in PyPSA (i.e. 230_AAL_DSS, 230=nominal_voltage, AAL=unique substation name in CODERS, DSS=Substation type)
    return bus_x: Bus longitude. 
    return bus_y: Bus latitude.
    '''
    # (1) identical match
    for idx,substation in df_substations.iterrows():
        if substation["node_code"] == node_code:
            bus_name = str(line["v_nom"]) + "_" + "_".join(node_code.split('_')[1:]) # (i.e. 230_AAL_DSS)
            bus_x = substation["longitude"]
            bus_y = substation["latitude"]
            return bus_name, bus_x, bus_y

    # print(f"Did not find exact match for line node: {node_code}") # To-be logged

    # (2) International and Interprovincial nodes
    if node_code.split('_')[-1] in INTERTIE_SUFFIXES:
        bus_name = str(line["v_nom"]) + "_" + "_".join(node_code.split('_')[1:])
        if node_code == "PP_BCAB3_IPT":
            # ~ 50 km east
            bus_y, bus_x = 50.247937, -114.2 
        elif node_code == "PP_BCAB1_IPT":
            # ~ 18 km east
            bus_y, bus_x = 49.735535, -114.6
        elif node_code == "PP_BCAB4_IPT":
            # ~21 km east
            bus_y, bus_x = 58.64525, -119.7
        elif node_code == "XX_BCUS2_INT":
            # ~ 1 km south
            bus_y, bus_x = 48.9974, -117.341514
        elif node_code == "XX_BCUS1_INT":
            # ~ 21 km south
            bus_y, bus_x = 48.97 , -122.873948
        elif node_code == "PP_BCAB2_IPT":
            # ~ 108 km east
            bus_y, bus_x = 49.500543, -114.08
        elif node_code == "PP_ABSK1_IPT":
            #
            bus_y, bus_x = 51.283272, -107.427622
        elif node_code == "XX_ABUS1_INT":
            #
            bus_y, bus_x = 48.312902, -112.660219
        elif node_code == "XX_SKUS1_INT":
            bus_y, bus_x = 48.720572, -105.863443

        elif node_code == "PP_MBSK1_IPT":
            bus_y, bus_x = 54.794729, -101.885937

        elif node_code == "PP_MBSK2_IPT":
            bus_y, bus_x = 53.722668, -101.769565

        elif node_code == "PP_MBSK3_IPT":
            bus_y, bus_x = 49.495283, -101.393102

        elif node_code == "PP_MBSK4_IPT":
            bus_y, bus_x = 50.544463, -101.474715
        
        elif node_code == "PP_MBSK5_IPT":
            bus_y, bus_x = 51.202531, -101.538936

        else:
            log.error('no match for IPT/INT substation {}'.format(node_code))
            exit(3)
        return bus_name, bus_x, bus_y

    # (3) Look in the generator dataset for a location for the node_code
    # First match used for location
    for _,generator in generators.iterrows():
        if generator["connecting_node_code"] == node_code:
            bus_name = str(line["v_nom"]) + "_" + "_".join(node_code.split('_')[1:]) # (i.e. 230_AAL_DSS)
            bus_x = generator["longitude"]
            bus_y = generator["latitude"]
            return bus_name, bus_x, bus_y
        
    # (4) Find first matching 3-middle characters
    for idx,substation in df_substations.iterrows():
        if node_code.split('_')[1] == substation['node_code'].split('_')[1]:
            bus_name = str(line["v_nom"]) + "_" + "_".join(node_code.split('_')[1:])
            bus_x = substation["longitude"]
            bus_y = substation["latitude"]
            return bus_name, bus_x, bus_y
    # print(f"Did not find partial match for: {node_code}") # To-be logged


    # (6) Find first matching 3-middle characters
    # First match used for location
    for _,generator in generators.iterrows():
        if node_code.split('_')[1] == generator["connecting_node_code"].split('_')[1]:
            bus_name = str(line["v_nom"]) + "_" + "_".join(node_code.split('_')[1:]) # (i.e. 230_AAL_DSS)
            bus_x = generator["longitude"]
            bus_y = generator["latitude"]
            return bus_name, bus_x, bus_y

    # (5) No match found...
    # Case of no matches found....
    log.warning(f"no information to create bus for: {node_code}")

    return None,None,None

def sanitize(df_substations:pd.DataFrame)->pd.DataFrame:
    """
    This function sanitizes the input DataFrame by ensuring that the 'summer_rating_mw' column is numeric 
    and that the 'starting_node_code' and 'ending_node_code' columns are stripped of whitespace.
    
    Parameters:
    df_substations (pd.DataFrame): The input DataFrame containing substation data.
    
    Returns:
    pd.DataFrame: The sanitized DataFrame with cleaned 'summer_rating_mw', 'starting_node_code', and 'ending_node_code' columns.

    Raises:
        AssertionError: If the input is not a DataFrame or if required columns are missing.
    """
    assert isinstance(df_substations, pd.DataFrame), "Input must be a pandas DataFrame"
    assert "node_code" in df_substations.columns, "DataFrame must contain 'node_code' column"

    # Strip whitespace from 'starting_node_code' and 'ending_node_code' columns
    node_cols = REQUIRED_COLUMNS
    for col in node_cols:
        assert col in df_substations.columns, f"DataFrame must contain '{col}' column"
        # Option A: Apply .str accessor directly (pandas preserves NaNs natively)
        df_substations[col] = df_substations[col].str.replace(r"\s+", "", regex=True)
        
    

    return df_substations

def get_prepared_buses(prepared_lines:pd.DataFrame=None)->pd.DataFrame:
    '''
    This script prepares the csv files for creating the base PyPSA_BC network.
    # Outfiles: 
    - buses.csv, 
    - lines.csv, 
    - line_types.csv, 
    - transformers.csv, 
    - transformer_types.csv

    '''
    utils.print_update(level=1,message="Preparing base nework for PyPSA_BC")

    # A) load data
    df_sub = coders_data.substations
    prepared_lines_df:pd.DataFrame = prepared_lines if prepared_lines is not None else lines.get_prepared_lines()[0]
    generators_df = coders_data.generators

    # B) process data
    df_sub = sanitize(df_sub)

    # # (5) Create dataframe of BC buses from the lines and substations
    df_buses_bc = create_bus_df(prepared_lines=prepared_lines_df, 
                                prepared_substations=df_sub, 
                                prepared_generators=generators_df)
    
    

    """ 
    #--------------------------------------------------------------------------------
    skipping, EL_20260723
    
    # ######### MODIFICATIONS: To adjust for the CODERS updates in naming. #####################
    # lines + buses modified (NOTE: lines name is NOT modified)
    # i) add all other mismatches between hydro,wind,sol,tpp generators and the buses.csv file
    # Find all "_DFS" endings and correct for each.
    hydro_codes = {'230_CMS_DFS':'230_CMS_GSS', '138_JOR_DFS':'138_JOR_GSS','138_PUN_DFS':'138_PUN_GSS',
                   '69_LB1_DFS':'69_LB1_GSS', '69_SON_DFS':'69_SON_GSS','69_LAJ_DFS':'69_LAJ_GSS',
                   '69_SPN_DFS':'69_SPN_GSS','230_RGA_DSS':'230_RGA_TSS'} # CRS and SPN not found.
    
    utils.print_update(level=2,message="Checking and calibrating hydro generator data...")
    utils.fix_coders_update(df_buses_bc, col_to_correct='name',codes=hydro_codes)
    utils.fix_coders_update(lines, col_to_correct='bus0',codes=hydro_codes)
    utils.fix_coders_update(lines, col_to_correct='bus1',codes=hydro_codes)

    utils.print_update(level=2,message="Checking and calibrating bus data...")
    check_missing_buses(df_sub, lines)
    #--------------------------------------------------------------------------------
     
    """
    check_missing_buses(prepared_substations=df_sub, prepared_lines=prepared_lines_df)
    # # Additional attributes
    # df_buses_bc['substation_type'] = df_buses_bc['name'].apply(lambda x: x.split('_')[-1])
    buses = df_buses_bc 
    
    buses.to_csv(Path(aparser.data_cfg.get('output').get('base_network'))/"buses.csv", index=False)
    log.debug(f"buses.csv saved to {Path(aparser.data_cfg.get('output').get('base_network'))/'buses.csv'}")
    return buses