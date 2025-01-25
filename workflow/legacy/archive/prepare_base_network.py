import pandas as pd
from bc_combined_modelling import utils
import math


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

def check_missing_buses(df_sub_bc, df_lines_bc):
    '''
    Checks for substations which are missing from the lines dataset.
    The only bus like this for BC is BC_WAX_GSS which should be connceted to BC_WAN.
    df_buses_bc: DF being prepared for saving and loading into PyPSA.
    df_sub_bc: DF from CODERS of all substations.
    '''
    # Checking to BC_WAX_GSS is the only one missing
    sub_unique_codes = df_sub_bc["node_code"].apply(lambda x: x.split('_')[1].strip(' ')).unique().tolist()
    line_unique_codes_start = df_lines_bc["starting_node_code"].apply(lambda x: x.split('_')[1].strip(' ')).unique().tolist()
    line_unique_codes_end = df_lines_bc["ending_node_code"].apply(lambda x: x.split('_')[1].strip(' ')).unique().tolist()
    line_unique_codes = list(set(line_unique_codes_start + line_unique_codes_end))

    for code in sub_unique_codes:
        if code not in line_unique_codes:
            print(f"the code {code} is missing from the lines dataset!") # to be logged

def add_pypsa_columns_2_line_df(df):
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
        voltage_type = f'{str(int(line["voltage_in_kv"]))}kV'
        name.append(line["starting_node_code"][3:] + line["ending_node_code"][2:])
        type.append(voltage_type)
        bus0.append(voltage_type.rstrip('kV') + "_" + line['starting_node_code'].split('_')[1] + "_" + line['starting_node_code'].split('_')[2])
        bus1.append(voltage_type.rstrip('kV') + "_" + line['ending_node_code'].split('_')[1] + "_" + line['ending_node_code'].split('_')[2])
        length.append(line["line_segment_length_km"])
        v_nom.append(line["voltage_in_kv"])

    df['name'] = name
    df['type'] = type
    df['bus0'] = bus0
    df['bus1'] = bus1
    df['length'] = length
    df['v_nom'] = v_nom

def create_bus_df(df_lines, df_substations, generators):
    '''
    This function will create an initial DataFrame of buses for PyPSA_BC from a DataFrame of lines.
    When creating the buses it
    The lines DF contains the node names for the buses, nomial voltage, and the carrier is implicitly added.
    '''
    # name = []
    # x = []
    # y = []
    # type = []
    # v_nom = []
    data_dict = {}

    # (1) Add buses based on line nodes
    for _,line in df_lines.iterrows():
        # Search for match between line and substation
        for node_code in [line["starting_node_code"], line["ending_node_code"]]:
            bus_name, bus_x, bus_y = get_bus_name_x_y(line, node_code, df_substations, generators)
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
    if node_code.split('_')[-1] in ["IPT","INT"]:
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
            print('ERROR: Did not find a match for IPT/INT substation {}'.format(node_code))
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
    print(f"There is no information to create bus for: {node_code}") # To-be logged

    return None,None,None

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
    f_nom = 60 # nominal frequency in NA is 60 Hz
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
    voltage_2_react = {63:0.3800, 66:0.3821, 69: 0.3843, 72:0.3864, 115: 0.4171, 120:0.4207,
                       138:0.4336, 144:0.4379, 161:0.4500, 230:0.4800, 240:0.4769, 287:0.4621, 315:0.4532,
                       345:0.4438, 360:0.4391, 500:0.3950, 735:0.3800, 765:0.3781}

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
                    print('ERROR: TRANSMISSION LINE MISSING A MAPPING FROM VOLTAGE TO REACTANCE')
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

def add_line_op_params(df_lines_bc, df_line_types_bc):
    '''
    This function add operational parameters to the lines such as:
    s_nom = which is pulled from 
    '''
    df_lines_bc['s_nom'] = df_lines_bc.apply(lambda line: add_line_s_nom(line),axis=1)

def add_line_s_nom(line):
    '''
    This function is applied row-wise to calculate the s_nom (MVA) for each line.
    Using data from CODERS to impute the s_nom value for lines as follows:
    rule 1: use summer_rating_in_mva
    rule 2: use summer_capacity_in_mw adjusted by a power factor = 0.9
    Using the summer values is a pessimistic assumption for other seasons.
    '''
    if line['summer_rating_in_mva'] == 0:
        s_nom = round(line['summer_capacity_in_mw'] / 0.9, 4) # (assumed 0.9 power factor + rounded to 4th decimal)
    else:
        s_nom = line['summer_rating_in_mva']

    return s_nom


def create_transformer_df(df_buses):
    '''
    This function creates a dataframe of transformers based on buses with multiple voltages.
    Voltages attached incrementally from low to highest at a given bus.
    Transformers are each given as a standardized type.
    df_buses: Dataframe of PyPSA formatted buses.
    '''
    # Find unique buses and their indices
    bus_dict = {}
    for idx,row in df_buses.iterrows():
        key = "_".join(row['name'].split('_')[1:]) # i.e. VIT_TSS
        if key not in bus_dict.keys():
            bus_dict[key] = [row['v_nom']]
        else:
            bus_dict[key].append(row['v_nom'])

    #
    transformers = [] # list to hold transformers to create
    for bus,voltages in bus_dict.items():
        if len(voltages) <= 1: 
            continue
        else: # More than 1 voltage at the unique bus
            N = len(voltages)
            voltages_sorted = sorted(voltages)
            for idx in range(N-1):
                hv = voltages_sorted[idx+1]
                lv = voltages_sorted[idx]
                bus0 =  str(hv) + "_" + bus
                bus1 =  str(lv) + "_" + bus
                type = f"{hv}/{lv}"
                transformer_name = f'{bus}_{hv}_{lv}'
                transformers.append([transformer_name,
                                    bus0,
                                    bus1,
                                    type]
                                    )

    df_transformers = pd.DataFrame(transformers, columns = ['name','bus0','bus1','type'])

    return df_transformers
    
def create_tranformer_types_df(df_transformers):
    '''
    This function will create transformers for typical hv to lv lines.
    Assumption 1: All buses use standardized transformers and 
    there is only 1 type for each unique tuple of high and low voltage.
    Assumption 2: The capacity (MVA) of the transformer is not a bottleneck of the system,
    therefore, the capacity (s_nom) is set to be 2000 MVA (assumed limitless).
    Assumption 3: All transformers have the same parameters. (UPDATE LATER) 
    '''
    data = []
    for idx,row in df_transformers.iterrows():
        hv = int(row['type'].split('/')[0])
        lv = int(row['type'].split('/')[1])
        name = row['type']
        f_nom = 60
        s_nom = 2000
        v_nom_0 = hv
        v_nom_1 = lv
        vsc = 10 # Update
        vscr = 0.3 # Update
        pfe= 30 # Update
        i0 = 0.04 #2-10%
        phase_shift = 150 
        tap_side = 0
        tap_neutral = 0
        tap_min = -9
        tap_max = 9
        tap_step = 1.5

        data.append([name,f_nom,s_nom,v_nom_0,v_nom_1,
                    vsc,vscr,pfe,i0,phase_shift,tap_side,
                    tap_neutral,tap_min,tap_max,tap_step])

    df_transformer_types = pd.DataFrame(data,columns=["name","f_nom","s_nom",
                               "v_nom_0","v_nom_1", "vsc",
                               "vscr","pfe","i0","phase_shift",
                               "tap_side","tap_neutral","tap_min",
                               "tap_max","tap_step"]).drop_duplicates()

    return df_transformer_types



def main():
    '''
    This script prepares the csv files for creating the base PyPSA_BC network.
    outfiles: buses.csv, lines.csv,, line_types.csv, transformers.csv, transformer_types.csv

    '''

    # Read in configuration file
    config_file = r"config/pypsa_config.yaml"
    cfg = utils.load_config(config_file)

    # A) load data
    # /mnt/c/Users/pmcw9/Delta-E/PICS/Data
    transmission_line_path = cfg['data']["coders"]["lines"]
    substations_path = cfg['data']["coders"]["substations"]
    generators_path = cfg['data']['coders']['generators']
    transmission_line_type_table = cfg['data']["custom"]["line_table"] # NOTE: To be removed eventually
    

    df_lines = pd.read_csv(transmission_line_path)
    df_substations = pd.read_csv(substations_path)
    generators = pd.read_csv(generators_path)
    df_line_table = pd.read_excel(transmission_line_type_table) # Tables with line type data for indexing by ampacities

    # Loop over all regions to include
    df_sub = df_substations[df_substations["province"].apply(lambda x: x in cfg['output']['prepare_base_network']['regions'])].copy()
    df_lines = df_lines[df_lines["province"].apply(lambda x: x in cfg['output']['prepare_base_network']['regions'])].copy()
    
    # MODIFICATIONS 2024-10-01: Simialr to create_hydro_asset fix for the bridge cascade. Here the substation name for the following is modified:
    # BC_BR1_DFS modified to BC_BR1_GSS
    bridge_codes = {'BC_BR1_DFS':'BC_BR1_GSS'}
     
    # i)
    mask = df_sub['node_code'].isin(bridge_codes)
    df_sub.loc[mask,'node_code'] = df_sub.loc[mask,'node_code'].apply(lambda old_cold: bridge_codes[old_cold])
    # ii)
    mask = df_lines['starting_node_code'].isin(bridge_codes)
    df_lines.loc[mask,'starting_node_code'] = df_lines.loc[mask,'starting_node_code'].apply(lambda old_cold: bridge_codes[old_cold])
    # iii) 
    mask = df_lines['ending_node_code'].isin(bridge_codes)
    df_lines.loc[mask,'ending_node_code'] = df_lines.loc[mask,'ending_node_code'].apply(lambda old_cold: bridge_codes[old_cold])

    # B) process data

    # (0) Replace NaN for summer rating with 0
    df_lines['summer_rating_in_mva'] = df_lines["summer_rating_in_mva"].fillna(0.)

    # (1) Correction to data
    correct_line_node_name(df_lines) 

    # (2) Remove spaces from code names
    df_sub["node_code"] = df_sub["node_code"].apply(lambda x: x.replace(" ",""))
    df_lines["starting_node_code"] = df_lines["starting_node_code"].apply(lambda x: x.replace(" ",""))
    df_lines["ending_node_code"] = df_lines["ending_node_code"].apply(lambda x: x.replace(" ",""))

    # (3) Add missing lines to dataset
    lines = df_lines
    # NOTE: For BC only
    lines = add_missing_lines(df_lines)

    # (4) Enrich coders dataframe of BC lines with columns used by PyPSA 
    add_pypsa_columns_2_line_df(lines)

    # # (5) Create dataframe of BC buses from the lines and substations
    df_buses_bc = create_bus_df(lines, df_sub, generators) 
    
    

    # ######### MODIFICATIONS: To adjust for the CODERS updates in naming. #####################
    # lines + buses modified (NOTE: lines name is NOT modified)
    # i) add all other mismatches between hydro,wind,sol,tpp generators and the buses.csv file
    # Find all "_DFS" endings and correct for each.
    hydro_codes = {'230_CMS_DFS':'230_CMS_GSS', '138_JOR_DFS':'138_JOR_GSS','138_PUN_DFS':'138_PUN_GSS',
                   '69_LB1_DFS':'69_LB1_GSS', '69_SON_DFS':'69_SON_GSS','69_LAJ_DFS':'69_LAJ_GSS',
                   '69_SPN_DFS':'69_SPN_GSS','230_RGA_DSS':'230_RGA_TSS'} # CRS and SPN unfound.
    utils.fix_coders_update(df_buses_bc, col_to_correct='name',codes=hydro_codes)
    utils.fix_coders_update(lines, col_to_correct='bus0',codes=hydro_codes)
    utils.fix_coders_update(lines, col_to_correct='bus1',codes=hydro_codes)



    check_missing_buses(df_sub, lines)

    # # (6) create dataframe of line types for BC
    bc_line_types = create_line_types_df(lines, df_line_table)

    # # (7) rename lines which are duplicates (add suffix of _#)
    rename_duplicate_lines(lines)

    # # (8) add all needed operational parameters to lines
    add_line_op_params(lines, bc_line_types)

    # # (9) create dataframe of transformers for BC
    df_transformers_bc = create_transformer_df(df_buses_bc)

    # # (10) create dataframe of transformer types for BC
    df_transformer_types_bc = create_tranformer_types_df(df_transformers_bc)

    # # Additional attributes
    df_buses_bc['substation_type'] = df_buses_bc['name'].apply(lambda x: x.split('_')[-1])
    buses = df_buses_bc 

    # C) record data
    path = cfg['output']['prepare_base_network']['folder']
    utils.create_folder(path)
    lines.to_csv(path + "/lines.csv", index=False,
                    columns=['name','type','bus0','bus1','length','v_nom','s_nom'])
    buses.to_csv(path + "/buses.csv", index=False)
    bc_line_types.to_csv(path + "/line_types.csv", index=False)
    df_transformers_bc.to_csv(path + "/transformers.csv", index=False)
    df_transformer_types_bc.to_csv(path + "/transformer_types.csv", index=False)


if __name__ == '__main__':
    main()