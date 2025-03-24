import pandas as pd
from pypsa_bc import utils
import warnings
from pathlib import Path
from typing import Optional
import datetime
from models.PyPSA_BC.workflow.scripts import (
    prepare_base_network,
    create_hydro_assets,
    create_reservoir_inflows,
    create_ror_ps,
    create_ext_wind_assets,
    create_ext_wind_ts,
    create_ext_solar_assets,
    create_ext_solar_ts,
    create_ext_tpp_assets,
    create_cutout,
    enrich_format_hydro,
    enrich_format_tpp,
    enrich_format_vre,
    disaggregate_load,
    build_model
)

# Suppress specific warnings
warnings.filterwarnings("ignore", category=FutureWarning)
""" 
pypsa_model_build_args={
    'copperplate':False,
    'update_data':False,
    'update_load':True,
    'year':2021,
    'resource_options': 'investment',  # 'investment' or 'full_potential',
    'solved_network_save_to': Path('results/pypsa')
}
"""
# -------------------------------------------------------------------------------------------------------------------------------------
def run_data_preparation_workflow():
    prepare_base_network.main()
    create_hydro_assets.main()
    create_reservoir_inflows.main()
    create_ror_ps.main()
    # create_cutout.main() # time intensive job
    create_ext_wind_assets.main()
    create_ext_wind_ts.main()
    create_ext_solar_assets.main()
    create_ext_solar_ts.main()
    create_ext_tpp_assets.main()

def enrich_format_workflow():
    enrich_format_hydro.main()
    enrich_format_tpp.main()
    enrich_format_vre.main()

def main(copperplate:bool, # the EV load data is prepared for Copperplate (single region)
        update_data:bool,
        update_load:bool,
        resource_options:str,
        start_date:str,
        end_date:str,
        ev_charging:str,
        ev_population:float,
        run_tag:Optional[int]=None,
        solved_network_save_to:str= 'results/pypsa'):
    """ 
    Args:
        update_data (bool) : 'True ' or 'False
        update_load
        copperplate (bool) : 'True ' or 'False
        'resource_options' (str): 'investment' or 'full_potential'
        year (int): 2021 to 2050
    """
    year:int= datetime.datetime.strptime(start_date, "%Y-%m-%d").year
    
    if update_data:
        utils.print_update(level=1, message='Initiating network and profile data preparation...')
        run_data_preparation_workflow()    
        enrich_format_workflow()

    elif update_load:
        utils.print_update(level=1, message='Skipping network and profile data update. Using the prepared data.')
        
    else:
        utils.print_update(level=1, message='Skipping network and profile data update. Using the prepared data from <data/pypsa/pypsa_data>')
        utils.print_update(level=1, message='Skipping load data update. Using the prepared data.')

    load_bch_raw:pd.DataFrame=pd.read_excel('data/pypsa/downloaded_data/load/bc_hydro_load/BalancingAuthorityLoad2021.xls')
    provincial_load_MWh:float=disaggregate_load.fix_hourly_load(load_bch_raw,year)

    # replace with any total load data
    provincial_total_load_MWh:float=provincial_load_MWh.LOAD.sum() #MWh, provincial total
    disaggregate_load.main(provincial_total_load_MWh)
    
    build_main_args= {
    'copperplate':copperplate,
    'capacity_choice':resource_options,
    'year':year,
    'ev_charging': ev_charging,
    'ev_population':ev_population,
    'run_tag':run_tag,
    }
    build_model.main(**build_main_args)


if __name__ == '__main__':
    main()