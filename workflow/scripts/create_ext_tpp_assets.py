from pypsa_bc import utils
import pandas as pd

# handles the config loading centrally
from pypsa_bc.attributes_parser import AttributesParser
pypsa_aparser=AttributesParser()


def get_fuel_bus(site, fuel_type, cfg):
    '''
    Returns name of the fuel bus for a specific tpp.
    '''
    if cfg["output"]["create_ext_tpp_assets"]["gas_grid"]:
        fuel_bus = "{} {} Bus".format(site.name, fuel_type)
    else:
        fuel_bus = "Global {} Bus".format(fuel_type)

    return fuel_bus

def get_tpp_dict(site, gen_params, bus_dict, tpp_gen_types, cfg, province ='Not Selected'):
    '''
    Creates a thermal power plant dictionary that is used to add network components of the thermal Powerplants.
    NOTE: Assume running for hourly (need to update later)
    NOTE: Fuel assumed to be in units of MMBtu
    '''
    fuel_type = tpp_gen_types[site["gen_type"]]

    fuel_bus = get_fuel_bus(site, fuel_type, cfg)
    elc_bus = utils.get_gen_bus(site['connecting_node_code'], bus_dict, province)

    name = " ".join([site.name, site["gen_type"], "Link"])

    # Efficiency inverse of the heat rate
    eff_hr = 1. / gen_params["heat_rate"]

    # fuel_cost + variable_cost_MWH * conversion in terms of MMBtu (var_cost * efficiency + fuel_cost)
    # NOTE: Cost for the fuel could be listed here or assigned to production method on the NG bus. 
    marginal_cost = gen_params["average_fuel_price_CAD_per_MMBtu"] + gen_params["variable_om_cost_CAD_per_MWh"] * eff_hr
    
    # Create link + store representation of generator
    # bus_name = " ".join([fuel_type, "Bus"])
    tpp_comp_dict = {"class_name":"Link",
                    "name": name,
                    "bus0": fuel_bus,
                    "bus1": elc_bus,
                    "carrier": fuel_type,
                    "efficiency":eff_hr,
                    "ramp_limit_up":min(gen_params["ramp_rate_percent_per_min"]*60, 1), #* site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "ramp_limit_down":min(gen_params["ramp_rate_percent_per_min"]*60, 1), #* site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "p_nom_extendable":False,
                    "committable":cfg["output"]["create_ext_tpp_assets"]["UC"],
                    "min_up_time":gen_params["min_up_time_hours"],
                    "min_down_time":gen_params["min_down_time_hours"],
                    # "ramp_limit_start_up":row["ramp_limit_start_up"], # no data atm
                    # "ramp_limit_shut_down":row["ramp_limit_shut_down"], # no data atm
                    "p_nom":site["facility_installed_capacity"] * gen_params["heat_rate"],
                    "marginal_cost":marginal_cost,
                    # "p_min_pu":gen_params["min_plant_load"] # watch out for the forced run condition when UC is off.
                    }

    return tpp_comp_dict

def write_tpp_csv(tpp_assets, gen_generic, bus_dict, tpp_gen_types, cfg, province='Not Selected'):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing TPP facilities in PyPSA.
    '''
    tpp_dict = {}

    # Next create the dictionaries for PyPSA instantiation
    for _,site in tpp_assets.iterrows():

        fuel_type = tpp_gen_types[site["gen_type"]]

        # if bool(cfg["output"]["enrich_format_tpp"]["gas_grid"]): # IF THERE IS A GAS GRID ONLY
        #     fuel_bus = "{} {} Bus".format(site['asset_id'], fuel_type)

        # else: # GLobal gas store will be used...
        #     fuel_bus = "Global {} Bus".format(fuel_type)

        #     if fuel_bus not in tpp_dict.keys():
        #         tpp_dict['Global NG Bus'] = {"class_name": "Bus",
        #                                      "name":fuel_bus,
        #                                      "carrier":fuel_type}
        #     else:
        #         pass

    # Add fuel_bus
        aid = site.name
        gen_params = gen_generic[gen_generic["generation_type"] == site["gen_type"]].squeeze()
        tpp_dict[aid] = get_tpp_dict(site, gen_params, bus_dict, tpp_gen_types, cfg, province)

    # write csv
    out_file = cfg["output"]["create_ext_tpp_assets"]["fname"]
    df_tpp_assets = pd.DataFrame.from_dict(tpp_dict, orient='index').reset_index()
    
    df_tpp_assets.to_csv(out_file, index=False)

def main():
    '''
    This script creates the dictionaries needed to instantiate the thermal power plants (TPP) in PyPSA_BC.
    '''
    # get configuration file
    cfg=pypsa_aparser.pypsa_cfg

    print(f"Preparing existing thermal power assets...")


    gen_generic = pd.read_csv(cfg["data"]["coders"]["gen_generic"])
    gens = pd.read_csv(cfg["data"]["coders"]["generators"])
    # hist_gen = pd.read_csv(cfg["output"]["create_tpp_assets"]["cogen_history"],parse_dates=True, index_col=0)
    buses = pd.read_csv(cfg['output']['prepare_base_network']['folder'] + "/buses.csv")['name'].tolist()

    # All generation types which are thermal PP in the CODERS dataset.
    tpp_gen_types = {'NG_CT':"NG", 'NG_CC':"NG", 'gasoline_CT':"NG",
                    'oil_CT':"Oil", "coal":"Coal", 'oil_ST':"Oil",
                    'diesel_CT':"Diesel", 'coal_CCS':"Coal", 'biomass':'Biomass',
                    "biogas":"Biomass", "Coal_CCS":"Coal",
                    'NG_CG':'NG'}

    # Determine tpp generators in BC
    
    # In BC all generators are "NG" however in general this is not true.
    # NOTE: For time being there will only be a single NG bus. however, in the future.
    # Buses will need to be added for each node.


    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    bus_dict = utils.create_standard_gen_bus_map(buses)
    
    # (1) Write pickle dictionaries for the tpp assets.
    # per fuel type now to access specific costs for each

    mask = (gens['province'] == 'BC') & (gens["gen_type"].apply(lambda x: x in tpp_gen_types ))
    subset = ["connecting_node_code"]
    sum_list = ["facility_installed_capacity","facility_average_annual_energy"]
    gens[mask].groupby(subset, group_keys=False).apply(lambda x: utils.merge_assets(x, subset, sum_list))


    for region in cfg['output']['prepare_base_network']['regions']:
        mask = (gens["province"] == region) & (gens["gen_type"].apply(lambda x: x in tpp_gen_types ))
        tpp_gens = gens[mask].copy()

        # Code to aggregate tpp units
        subset = ["connecting_node_code"]
        sum_list = ["facility_installed_capacity","facility_average_annual_energy"]
        tpp_gens['node_code'] = tpp_gens['connecting_node_code']
        tpp_gens = tpp_gens[mask].groupby(subset, group_keys=False).apply(lambda x:
                                        utils.merge_assets(x, subset, sum_list))
        # Code to aggregate tpp units
        tpp_gens["connecting_node_code"] = tpp_gens.index
        tpp_gens.set_index('node_code', inplace=True)

        ext_tpp_assets = utils.add_generic_columns_tpp(tpp_gens, gen_generic)

        # write_tpp_csv(tpp_gens, gen_generic, bus_dict, tpp_gen_types, cfg, province = region)

    # Modified 2024-10-04: Fix for the node codes with no matching buses
    codes = {"BC_CRS_DSS":"BC_CRS_DFS", "BC_DGB_TSS":"BC_DGB_DSS"}
    utils.fix_coders_update(ext_tpp_assets,'connecting_node_code',codes)
    ext_tpp_assets.to_csv(cfg['output']['create_ext_tpp_assets']['fname'],index=False)
    utils.print_update(level=2,message="Finished preparing existing thermal power assets.")
    
if __name__ == '__main__':
    main()