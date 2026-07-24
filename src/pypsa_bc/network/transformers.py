import math
import pandas as pd
from pypsa_bc import utils
from pypsa_bc.attributes_parser import AttributesParser
from pypsa_bc.data.coders import get_coders
from pypsa_bc.reporting.assumptions import log_assumption
from pypsa_bc.reporting.logger import get_logger
from pathlib import Path
from pypsa_bc.network import lines
from pypsa_bc.network import buses

log = get_logger("transformers")

# instance of Attributes Parser to access the base_network config and other attributes
aparser=AttributesParser()
BASE_NETWORK_CFG=aparser.base_network_cfg

# shared CODERSData instance (one per run; tables cached in-memory)
coders_data=get_coders()
REQUIRED_COLUMNS = coders_data.required.get('substations', [])
META = {"module":  Path(__file__).stem , 
        "scenario": aparser.get_scenario,
        "path": aparser.assumption_report_save_to}


def create_transformer_df(prepared_buses:pd.DataFrame=None)->pd.DataFrame:
    '''
    This function creates a dataframe of transformers based on buses with multiple voltages.
    Voltages attached incrementally from low to highest at a given bus.
    Transformers are each given as a standardized type.
    prepared_buses: Dataframe of PyPSA formatted buses.
    '''
    assert prepared_buses is not None, "prepared_buses cannot be None. Please provide a valid DataFrame of prepared buses."
    # Find unique buses and their indices
    bus_dict = {}
    skipped = []
    for idx,row in prepared_buses.iterrows():
        name = row['name']
        if not isinstance(name, str):          # NaN / non-string bus name
            skipped.append(idx)
            continue
        key = "_".join(name.split('_')[1:]) # e.g. "230_VIT_TSS" -> "VIT_TSS"
        if key not in bus_dict.keys():
            bus_dict[key] = [row['v_nom']]
        else:
            bus_dict[key].append(row['v_nom'])
    if skipped:
        log.warning(f"skipped {len(skipped)} bus row(s) with non-string name: {skipped}")
        log_assumption(
            parameter="transformer build — skipped nameless buses",
            value=f"{len(skipped)} bus row(s) without a valid name skipped",
            unit="-",
            rationale="Bus rows with NaN name (e.g. intertie boundary nodes) cannot form a transformer key",
            source="create_transformer_df",
            **META,
        )
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
    
def create_transformer_types_df(df_transformers):
    '''
    This function will create transformers for typical hv to lv lines.
    Assumption 1: All buses use standardized transformers and 
    there is only 1 type for each unique tuple of high and low voltage.
    Assumption 2: The capacity (MVA) of the transformer is not a bottleneck of the system,
    therefore, the capacity (s_nom) is set to be 2000 MVA (assumed limitless).
    Assumption 3: All transformers have the same parameters. (UPDATE LATER) 
    '''
    data = []
    transformer_attributes:dict = BASE_NETWORK_CFG.get('transformer',{})
    for idx,row in df_transformers.iterrows():
        hv = int(row['type'].split('/')[0])
        lv = int(row['type'].split('/')[1])
        
        name = row['type']
        f_nom = BASE_NETWORK_CFG.get('f_nom', 60) #hz
        s_nom = transformer_attributes.get('s_nom', 2000) # MVA
        v_nom_0 = hv 	# Nominal voltage, side 0 (HV)	kV	High-voltage winding.
        v_nom_1 = lv    # Nominal voltage, side 1 (LV)	kV	Low-voltage winding.
        vsc = transformer_attributes.get('vsc', 10) # %
        vscr = transformer_attributes.get('vscr', 0.3) # %
        pfe = transformer_attributes.get('pfe', 30) # kW
        i0 = transformer_attributes.get('i0', 0.04) # % '2-10% ; No-load / magnetising current ;	Sets shunt susceptance.'
        phase_shift = transformer_attributes.get('phase_shift', 150) # degrees
        tap_side = transformer_attributes.get('tap_side', 0) # 0/1; Which side carries the tap changer ; 0 = HV side.
        tap_neutral = transformer_attributes.get('tap_neutral', 0) # Tap position at nominal ratio; Neutral = no voltage adjustment.
        tap_min = transformer_attributes.get('tap_min', -9) # Tap range ; 19 positions
        tap_max = transformer_attributes.get('tap_max', 9) # Tap range ; 19 positions
        tap_step = transformer_attributes.get('tap_step', 1.5) # % ; Voltage change per tap ;	±9 × 1.5% ⇒ ±13.5% regulation range.

        data.append([name,f_nom,s_nom,v_nom_0,v_nom_1,
                    vsc,vscr,pfe,i0,phase_shift,tap_side,
                    tap_neutral,tap_min,tap_max,tap_step])
        
        log_assumption(parameter="transformer impedance (type params)", 
                       value="vsc=10%, vscr=0.3%, i0=0.04%, pfe=30kW",
                       unit="mixed",
                       rationale="Default nameplate values — refine per unit",
                       source="base network prep",
                       **META)

    df_transformer_types = pd.DataFrame(data,columns=["name","f_nom","s_nom",
                               "v_nom_0","v_nom_1", "vsc",
                               "vscr","pfe","i0","phase_shift",
                               "tap_side","tap_neutral","tap_min",
                               "tap_max","tap_step"]).drop_duplicates()

    return df_transformer_types

def get_prepared_transformers(prepared_buses:pd.DataFrame=None)->tuple[pd.DataFrame:,pd.DataFrame]:
    '''
    This script prepares the csv files for creating the base PyPSA_BC network.
    # Outfiles: 
    - buses.csv, 
    - lines.csv, 
    - line_types.csv, 
    - transformers.csv, 
    - transformer_types.csv

    '''
    assert prepared_buses is not None, "prepared_buses cannot be None. Please provide a valid DataFrame of prepared buses."
    
    df_transformers= create_transformer_df(prepared_buses=prepared_buses)
    df_transformer_types = create_transformer_types_df(df_transformers)
    
    df_transformers.to_csv(Path(aparser.data_cfg.get('output').get('base_network'))/"transformers.csv", index=False)
    log.debug(f"transformers.csv saved to {Path(aparser.data_cfg.get('output').get('base_network'))/'transformers.csv'}")
    
    
    df_transformer_types.to_csv(Path(aparser.data_cfg.get('output').get('base_network'))/"transformer_types.csv", index=False)
    log.debug(f"transformer_types.csv saved to {Path(aparser.data_cfg.get('output').get('base_network'))/'transformer_types.csv'}")
    
    return df_transformers, df_transformer_types