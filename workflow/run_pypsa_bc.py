import warnings

import pandas as pd
from workflow.scripts import (
    build_model,
    create_ext_solar_assets,
    create_ext_solar_ts,
    create_ext_tpp_assets,
    create_ext_wind_assets,
    create_ext_wind_ts,
    create_hydro_assets,
    create_reservoir_inflows,
    create_ror_ps,
    disaggregate_load,
    # create_cutout,
    enrich_format_hydro,
    enrich_format_tpp,
    enrich_format_vre,
    prepare_base_network,
)

from pypsa_bc import utils

# Suppress specific warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def main(update_data:bool=False):
    if update_data:
        utils.print_update(level=1, message='Initiating network and profile data preparation...')
        
        prepare_base_network.main()
        create_hydro_assets.main()
        create_reservoir_inflows.main()
        create_ror_ps.main()
        create_ext_wind_assets.main()
        create_ext_wind_ts.main()
        create_ext_solar_assets.main()
        create_ext_solar_ts.main()
        create_ext_tpp_assets.main()
        enrich_format_hydro.main()
        enrich_format_tpp.main()
        enrich_format_vre.main()
        
    else:
        utils.print_update(level=1, message='Skipping network and profile data preparation. Using the prepared data.')
    
    load_bch_raw:pd.DataFrame=pd.read_excel('data/downloaded_data/load/bc_hydro_load/BalancingAuthorityLoad2021.xls')
    provincial_total_load_MWh:float=disaggregate_load.fix_hourly_load(load_bch_raw,2021)

    # replace with any total load data
    provincial_total_load_MWh:float=provincial_total_load_MWh.LOAD.sum() #MWh, provincial total
    disaggregate_load.main(provincial_total_load_MWh)
    
    build_model.main()

    
if __name__ == '__main__':
    main()