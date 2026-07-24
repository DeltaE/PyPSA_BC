import math
import pandas as pd
from pypsa_bc import utils
from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.reporting.assumptions import log_assumption
from pypsa_bc.reporting.logger import get_logger
from pathlib import Path

log = get_logger("lines")

# instance of Attributes Parser to access the base_network config and other attributes
aparser=AttributesParser()
BASE_NETWORK_CFG=aparser.base_network_cfg

# shared CODERSData instance (one per run; tables cached in-memory)
coders_data=get_coders()
REQUIRED_COLUMNS = coders_data.required.get('transmission_lines', [])
META = {"module":  Path(__file__).stem , 
        "scenario": aparser.get_scenario,
        "path": aparser.assumption_report_save_to}


def sanitize(df_lines)->pd.DataFrame:
    """
    This function sanitizes the input DataFrame by ensuring that the 'summer_rating_mw' column is numeric 
    and that the 'starting_node_code' and 'ending_node_code' columns are stripped of whitespace.
    
    Parameters:
    df_lines (pd.DataFrame): The input DataFrame containing transmission line data.
    
    Returns:
    pd.DataFrame: The sanitized DataFrame with cleaned 'summer_rating_mw', 'starting_node_code', and 'ending_node_code' columns.

    Raises:
        AssertionError: If the input is not a DataFrame or if required columns are missing.
    """
    assert isinstance(df_lines, pd.DataFrame), "Input must be a pandas DataFrame"
    assert "summer_rating_mw" in df_lines.columns, "DataFrame must contain 'summer_rating_mw' column"
    
    # Convert 'summer_rating_mw' to numeric, coercing errors to NaN, then fill NaN with 0.0 and convert to float
    df_lines.loc[:, "summer_rating_mw"] = (
        pd.to_numeric(df_lines["summer_rating_mw"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )

    # Strip whitespace from 'starting_node_code' and 'ending_node_code' columns
    node_cols = ["starting_node_code", "ending_node_code"]
    for col in node_cols:
        assert col in df_lines.columns, f"DataFrame must contain '{col}' column"
        # Option A: Apply .str accessor directly (pandas preserves NaNs natively)
        df_lines[col] = df_lines[col].str.replace(r"\s+", "", regex=True)

    return df_lines

def add_line_op_params(lines: pd.DataFrame)->pd.DataFrame:
    '''
    This function add operational parameters to the lines such as:
    s_nom = which is pulled from CODERS summer_rating_in_mva or summer_rating_mw adjusted by a power factor defined in config.yaml with 'base_network: 'default_pf' key (default 0.9).
    
    arguments:
    lines: pd.DataFrame of the transmission lines in BC.
    line_types: pd.DataFrame of the transmission line types in BC.
    
    returns:
    lines: pd.DataFrame of the transmission lines in BC with operational parameters added.
    
    raises:
    AssertionError: If the input is not a DataFrame or if required columns are missing. 
    '''
    # Apply row-wise safely
    lines["s_nom"] = lines.apply(add_line_s_nom, axis=1)
    log.debug("s_nom calculated for all lines; assumption report updated")
    
    return lines

def get_line_table():
    '''
    This function retrieves the line table from CODERS and prepares it for use in a PyPSA network.
    
    It sanitizes the data, adds necessary columns, and calculates the 's_nom' value for each line.
    
    Returns:
    pd.DataFrame: A DataFrame containing the prepared transmission lines data for PyPSA.    
    '''
    transmission_line_type_table = aparser.data_cfg.get('inventory').get('line_table')
    df_line_table = pd.read_excel(transmission_line_type_table) # Tables with line type data for indexing by ampacities
    return df_line_table

def create_line_types_df(df_lines, df_line_table):
    '''
    This function adds line type information for each line.
    Assuming all lines can have their transmission inferred on the basis of ampacity alone...
    # eventually will need a calculator based on short, medium, or long and voltage.
    r: resistance per length (Ohm per km)
    x: resistance/reactance per length (Ohm per km)
    c: shunt capacitance per length (nF per km)
    i: Nominal current (kA)
    cc: Cross section (mm^2)
    1) match based on closest match in the table
    2) match based on average for similar lines???
    '''
    ampacity_sel_col = "winter_ampacity" # NOTE: Switched from winter to summer
    f_nom = BASE_NETWORK_CFG['f_nom'] # nominal frequency in NA is 60 Hz
    
    data_dict = {'name':[],
                'f_nom':[],
                'r_per_length':[],
                'x_per_length':[],
                'c_per_length':[],
                'i_nom':[],
                'mounting':[],
                'cross_section':[]}
    
    line_type_col = []
    amp_cap_2_idx = {amp_cap:idx for idx,amp_cap in enumerate(df_line_table["approx_current_capacity"])}

    # Manual dictionary of Voltage (kV) -> Reactance (ohm/km)
    voltage_2_react = BASE_NETWORK_CFG['voltage_to_reactance']

    for _,row in df_lines.iterrows():
        ampacity = row[ampacity_sel_col] # Later should look into making this based on a timeseries.
        if not math.isnan(ampacity):
            name = str(int(ampacity)) 
            line_type_col.append(name) # Add line type for matching within PyPSA.

            if name not in data_dict["name"]:
                idx = sorted([(abs(amp_cap-ampacity), idx) for amp_cap,idx in amp_cap_2_idx.items()])[0][-1] # Find row index for closest ampacity in the table.
                data_dict["name"].append(name) # More descriptive name later (ampacity for now).
                data_dict["f_nom"].append(f_nom) # Hz
                data_dict['r_per_length'].append(df_line_table["resistance_ac_25_deg"].iloc[idx] / 1000)
                # data_dict['x_per_length'].append(df_line_table["x_l"].iloc[idx])
                if row["voltage_in_kv"] in voltage_2_react.keys():
                    data_dict['x_per_length'].append(voltage_2_react[row["voltage_in_kv"]])
                else:
                    # find nearest voltage and linearly interpolate the value
                    log.error('transmission line missing a voltage->reactance mapping')
                    exit(2)
                data_dict['c_per_length'].append(8.85) # Assumed based on VI-PyPSA.. Needs updating..
                data_dict['i_nom'].append(df_line_table["resistance_ac_25_deg"].iloc[idx] / 1000)
                data_dict['mounting'].append("ol")
                data_dict['cross_section'].append(df_line_table["cross_section_mm2"].iloc[idx] )

            else:
                continue # already have this line type added
        else:
            # use mode ampacity for same voltage type 
            if row["voltage_in_kv"] == 63: # replace 63 kV since no ampacity on it
                voltage = 69 
            elif row["voltage_in_kv"] == 161: # replace 161 kV since no ampacity for it
                voltage = 138 
            elif row["voltage_in_kv"] == 72:
                voltage = 69
            else:
                voltage = row["voltage_in_kv"]

            # if len(df_lines[(df_lines["voltage_in_kv"] == voltage) & (~df_lines[ampacity_sel_col].isnull())][ampacity_sel_col]) == 0:
            #     ampacity = name # NOTE: Just using last ampacity that was valid...
            # else:
            ampacity = df_lines[(df_lines["voltage_in_kv"] == voltage) & (~df_lines[ampacity_sel_col].isnull())][ampacity_sel_col].mode()[0]
                
            line_type_col.append(str(int(ampacity)))

    df_lines["type"] = line_type_col
    df_line_types = pd.DataFrame(data_dict)

    return df_line_types

def add_pypsa_attributes(df):
    ''' 
    This function will add the columns to the line df which will be imported into a pypsa network.
    line names assigned according to standard of voltage (i.e. 230_AAL, 230 = 230kV and AAL = middle 3 char of node_code).
    Line parameters such as reactance and resistance are imputed based on line type.
    
    name: Name of the transmission line, formatted as the starting and ending node code appended together. (i.e. XXX_GSS_YYY_DSS)
    type: Voltage level of the transmission line (i.e. 230kV).
    bus0: Name of the starting bus.
    bus1: Name of the ending bus.
    v_nom: Nominal voltage level (i.e. 230).

    '''
    name = []
    type = []
    bus0 = []
    bus1 = []
    length = []
    v_nom = []

    for idx,line in df.iterrows():
        voltage_type = f'{str(int(line["voltage_kV"]))}kV'
        name.append(line["starting_node_code"][3:] + line["ending_node_code"][2:])
        type.append(voltage_type)
        bus0.append(voltage_type.rstrip('kV') + "_" + line['starting_node_code'].split('_')[1] + "_" + line['starting_node_code'].split('_')[2])
        bus1.append(voltage_type.rstrip('kV') + "_" + line['ending_node_code'].split('_')[1] + "_" + line['ending_node_code'].split('_')[2])
        length.append(line["line_segment_length_km"])
        v_nom.append(line["voltage_kV"])

    df['name'] = name
    df['type'] = type
    df['bus0'] = bus0
    df['bus1'] = bus1
    df['length'] = length
    df['v_nom'] = v_nom
    return df

def add_line_s_nom(line):
    '''
    This function is applied row-wise to calculate the s_nom (MVA) for each line. Using data from CODERS to impute the s_nom value for lines as follows:
    rule 1: use `summer_rating_in_mva`
    rule 2: use `summer_rating_mw` adjusted by a power factor defined in `config.yaml` with 'base_network: 'defaults.pf' key (default 0.9)
    
    Note:
        Using the summer values is a pessimistic assumption for other seasons.
    '''
    candidate_cols = ['summer_rating_in_mva', 'summer_rating_mw']

    # Assert that at least one candidate column is present in df_lines
    assert any(col in line.index for col in candidate_cols), (
        f"CODERS | table | transmission_lines | Expected at least one rating column "
        f"from {candidate_cols}, but found none in DataFrame columns: {line.index.tolist()}"
    )
    
    if "summer_rating_in_mva" in line.index and pd.notna(line["summer_rating_in_mva"]):
        s_nom = line['summer_rating_in_mva']
    else:
        s_nom = round(line['summer_rating_mw'] / BASE_NETWORK_CFG.get('pf', 0.9), 4) # (assumed 0.9 power factor + rounded to 4th decimal)

    
    log_assumption(parameter="line s_nom", 
                   value="summer_rating_mw (ttc_summer)", 
                   unit="MW", 
                   rationale=f"CODERS Total Transfer Capability (TTC) used as thermal proxy, converted to MVA using default power factor ({BASE_NETWORK_CFG.get('pf', 0.9)}) from base_network.yaml",
                   **META)
    
    return s_nom

def get_prepared_lines(df_lines:pd.DataFrame=None)->tuple[pd.DataFrame, pd.DataFrame]:
    """
    This function retrieves the transmission lines data from CODERS and prepares it for use in a PyPSA network.
    
    It sanitizes the data, adds necessary columns, and calculates the 's_nom' value for each line.
    
    Parameters:
    df_lines (pd.DataFrame, optional): A DataFrame containing transmission line data. If not provided, the function will load the data from CODERS.
    
    Returns:
    pd.DataFrame: A DataFrame containing the prepared transmission lines data for PyPSA.    
    """
    assert df_lines is None or isinstance(df_lines, pd.DataFrame), "Input must be a pandas DataFrame or None (default to load from CODERS)"

    if df_lines is None:
        df_lines=coders_data.transmission_lines
    
    # Make sure all required columns are present in the DataFrame
    
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df_lines.columns]
    assert not missing_cols, f" CODERS | table | transmission_lines | Required columns missing: {missing_cols}"

    # Sanitize Data
    df_lines=sanitize(df_lines)
    df_lines:pd.DataFrame = add_line_op_params(df_lines)
    
    lines:pd.DataFrame = add_pypsa_attributes(df_lines)
    line_types:pd.DataFrame = get_line_table()
    
    lines.to_csv(Path(aparser.data_cfg.get('output').get('base_network'))/"lines.csv", index=False)
    log.debug(f"lines.csv saved to {Path(aparser.data_cfg.get('output').get('base_network'))/'lines.csv'}")
    line_types.to_csv(Path(aparser.data_cfg.get('output').get('base_network'))/"line_types.csv", index=False)
    log.debug(f"line_types.csv saved to {Path(aparser.data_cfg.get('output').get('base_network'))/'line_types.csv'}")

    return lines,line_types

"""
#--------------------------------------------------------------------------------------------------------
# depricated functions; EL_20260723
#---------------------------------------------------------------------------------------------------------

def add_missing_lines(df_lines_bc):
    '''
    This function adds lines to the dataset which are missing.
    df_lines_bc: Dataframe of the transmission lines in BC.
    '''
    # (1) Add line connecting BC_WAX_GSS to BC_SEL_TSS
    mask = df_lines_bc["transmission_line_id"] == 14552 
    drop_cols = ['transmission_line_id', 'line_length_km', 'line_segment_length_km',
                'line_segment_length_mi', 'line_length_mi', 'starting_node_name',
                'starting_node_code','ending_node_name', 'ending_node_code']
    # data_dict = {k:v[0] for k,v in df_lines_bc[mask].drop(drop_cols,axis=1).to_dict(orient='list').items()}
    data_dict = df_lines_bc[mask].drop(drop_cols,axis=1).to_dict(orient='list')
    data_dict['transmission_line_id'] = [14552999]
    data_dict['line_length_km'] = [10]
    data_dict['line_segment_length_km'] = [10]
    data_dict['starting_node_code'] = ["BC_WAX_GSS"]
    data_dict['ending_node_code'] = ["BC_SEL_TSS"]

    return pd.concat([df_lines_bc, pd.DataFrame(data_dict)]).reset_index()

def correct_line_node_name(df_lines):
    ''' 
    This function corrects nodes in the lines dataset from CODERS which were improperly named.
    '''
    
    node_cols = ["starting_node_code", "ending_node_code"]
    old_2_new = {"BC_WHO_JCT":"BC_WAH_JCT",
                 "BC_GST_JCT":"BC_MCK_JCT",
                 "BC_VSY_JCT":"BC_NOR_DSS"}
    for col in node_cols:
        for old_node,new_node in old_2_new.items():
            df_lines.loc[df_lines[col] == old_node, col] = new_node
            
def rename_duplicate_lines(df_lines_bc):
    '''
    Renamed duplicate lines, adding _# to reach in order found.
    '''
    indices = df_lines_bc[df_lines_bc.duplicated(subset=["name"],keep=False)].index.tolist()
    name_dict = {} # keep record of modifications
    for index in indices:
        name = df_lines_bc.loc[index,"name"]
        if name not in name_dict.keys():      
            name_dict[name] = 1
        else:
            name_dict[name] += 1
        df_lines_bc.loc[index,"name"] = name + "_" + str(name_dict[name])

#-------------------------------------------------------------------------------------------------------------------
"""