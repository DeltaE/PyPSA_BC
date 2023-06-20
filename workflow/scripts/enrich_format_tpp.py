import pypsa
from bc_power import utils, hydro
import json
import pandas as pd



def get_tpp_dict(site, gen_params, bus_dict, tpp_gen_types, cfg):
    '''
    Creates a thermal power plant dictionary that is used to add network components of the thermal Powerplants.
    NOTE: Assume running for hourly (need to update later)
    NOTE: Fuel assumed to be in units of MMBtu
    '''
    fuel_type = tpp_gen_types[site["gen_type"]]
    
    if cfg['tpp']["gas_grid"]:
        fuel_bus = "{} {} Bus".format(site.name, fuel_type)
    else:
        fuel_bus = "Global {} Bus".format(fuel_type)

    elc_bus = utils.get_gen_bus(site.name, bus_dict)

    name = " ".join([site.name, site["gen_type"], "Link"])
    marginal_cost = gen_params["average_fuel_price_USD_per_MMBtu"]
    
    # Create link + store representation of generator
    # bus_name = " ".join([fuel_type, "Bus"])
    tpp_comp_dict = {"class_name":"Link",
                    "name": name,
                    "bus0": fuel_bus,
                    "bus1": elc_bus,
                    "carrier": fuel_type,
                    "efficiency":1. / gen_params["heat_rate_MMBtu_per_MWh"], # Output = (1 / (heat_rate)) * (fuel in mmBTU)
                    "ramp_limit_up":min(gen_params["ramp_rate_percent_per_min"]*60, 1) * site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "ramp_limit_down":min(gen_params["ramp_rate_percent_per_min"]*60, 1) * site["install_capacity_in_mw"], # Aggregated units needs adjustments
                    "p_nom_extendable":False,
                    "committable":cfg['tpp']["UC"],
                    "min_up_time":gen_params["min_up_time_hours"],
                    "min_down_time":gen_params["min_down_time_hours"],
                    # "ramp_limit_start_up":row["ramp_limit_start_up"], # no data atm
                    # "ramp_limit_shut_down":row["ramp_limit_shut_down"], # no data atm
                    "p_nom":site["install_capacity_in_mw"] * gen_params["heat_rate_MMBtu_per_MWh"],
                    "marginal_cost":marginal_cost, # cost per input unit (Need to be careful when combining fuel cost and variable cost)
                    "p_min_pu":gen_params["min_plant_load"] # watch out for the forced run condition when UC is off.
                    }

    return tpp_comp_dict

def write_tpp_dict(tpp_assets, gen_generic, bus_dict, tpp_gen_types, cfg):
    '''
    This function writes a dictionary containing the information needed to create the components for
    existing vre facilities in PyPSA.
    '''
    tpp_dict = {}

    # NEED CODE TO ADD NG BUSES.. Will need to think about this.
    # NOTE: Likley should be this eventually, since possible to have NG without tpp.
    

    # First aggregate the units
    subset = ["connecting_node_code"]
    sum_list = ["install_capacity_in_mw","annual_avg_energy_unit_in_gwh/y"] # Modify this later perhaps
    tpp_agg = tpp_assets.groupby("connecting_node_code",group_keys=False).apply(lambda x:
                                             utils.merge_assets(x, subset, sum_list))

    # Next create the dictionaries for PyPSA instantiation
    for _,site in tpp_agg.iterrows():
        fuel_type = tpp_gen_types[site["gen_type"]]
        if cfg['tpp']["gas_grid"]:
            fuel_bus = "{} {} Bus".format(site['asset_id'], fuel_type)
        else:
            fuel_bus = "Global {} Bus".format(fuel_type)
            if fuel_bus not in tpp_dict.keys():
                tpp_dict['Global NG Bus'] = {"class_name": "Bus",
                                             "name":fuel_bus,
                                             "carrier":fuel_type}
            else:
                pass

    # Add fuel_bus
        aid = site.name
        gen_params = gen_generic[gen_generic["generation_type"] == site["gen_type"]].squeeze()
        tpp_dict[aid] = get_tpp_dict(site, gen_params, bus_dict, tpp_gen_types, cfg)

    # write pickle
    out_file = cfg["pypsa_dict"]["components"] + cfg["pypsa_dict"]["tpp"]
    utils.write_pickle(tpp_dict, out_file)

def main():
    '''
    This script creates the dictionaries needed to instantiate the thermal power plants (TPP) in PyPSA_BC.
    '''
    # Read in configuration file
    config_file = r"/home/pmcwhannel/repos/PyPSA_BC/config/config.yaml"
    cfg = utils.load_config(config_file)

    gen_generic = pd.read_csv(cfg["coders"]["gen_generic"])
    gens = pd.read_csv(cfg["coders"]["generators"])

    buses = pd.read_csv(cfg["network"]["folder"] + "/buses.csv")['name'].tolist()

    # All generation types which are thermal PP in the CODERS dataset.
    tpp_gen_types = {'NG_CT':"NG", 'NG_CG':"NG", 'NG_CC':"NG",
                     'Gas_CT':"NG", 'Oil_CT':"Oil", "Coal":"Coal", 'Oil_ST':"Oil", 'Diesel_CT':"Diesel",
                    'Coal_CCS':"Coal"}

    # Determine tpp generators in BC
    mask = (gens["province"] == "BC") & (gens["gen_type"].apply(lambda x: x in tpp_gen_types ))
    tpp_gens = gens[mask].copy()

    # In BC all generators are "NG" however in general this is not true.
    # NOTE: For time being there will only be a single NG bus. however, in the future.
    # Buses will need to be added for each node.

    # (0A) Create folders if they have not been created already
    utils.create_folder(cfg["pypsa_dict"]['components'])

    # (0B) Get bus_dict for mapping node codes to PyPSA_BC ELC buses
    bus_dict = utils.create_standard_gen_bus_map(buses)

    # (0C) Get vre cost information
    
                                
    # (1) Write pickle dictionaries for the vre assets.
    write_tpp_dict(tpp_gens, gen_generic, bus_dict, tpp_gen_types, cfg)

    
if __name__ == '__main__':
    main()