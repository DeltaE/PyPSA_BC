import sys
from pathlib import Path
import pandas as pd


#Format the hourly load data to account for inconsistencies with how BC Hydro handles daylight savings
def fix_hourly_load(load, year):
    #Filter down to just the loads
    fixed = pd.to_numeric(load[load.columns[-1]], errors='coerce') #Rightmost column in BC Hydro hourly load spreadsheet, to_numeric converts the labels into NaN
    fixed = fixed.reset_index().drop(columns=['index'])
    fixed.columns = ['LOAD']

    #Remove NaN values (the column labels and sometimes the value for DST hour) and 0 values (sometimes the DST hour), then convert from MWh to kWh
    fixed = fixed.loc[fixed['LOAD'] > 0] * 1000
    
    #Index hourly loads by hourly timestamp
    fixed['TIME'] = pd.date_range(start=str(year)+'-01-01 00:00:00', end=str(year)+'-12-31 23:00:00', freq='h')
    fixed = fixed.set_index('TIME')

    #Return the hourly load for the entire province for a given year
    return fixed


#This part disaggregates the hourly load data from BC Hydro into 27 subdivisions of BC
def disaggregate(ceei, hourly):
    #Get all the regions of BC (order is retained)
    regions = ceei.loc[ceei['ORG_PART'] == 9000000]['ORG_NAME'].unique()

    #Total load for all of BC (Res + CSMI)
    ceei_total = ceei.loc[(ceei.ORG_NAME == 'British Columbia')]['CONSUMPTION_TOTAL'].sum()

    #Get proportion of annual loads per region
    proportions = pd.DataFrame(columns=['REGION', 'PROPORTION_RES', 'PROPORTION_CSMI'])

    #Create proportions dataframe
    for r in regions:
        #Skip the region for all of BC
        if r == 'British Columbia':
            break

        #Obtain proportion for residential and industrial separately
        proportion_res = ceei.loc[(ceei['ORG_NAME'] == r) & (ceei['ORG_PART'] == 9000000) & (ceei['SUB_SECTOR'] == 'Res')]['CONSUMPTION_TOTAL'].values[0] / ceei_total
        proportion_csmi = ceei.loc[(ceei['ORG_NAME'] == r) & (ceei['ORG_PART'] == 9000000) & (ceei['SUB_SECTOR'] == 'CSMI')]['CONSUMPTION_TOTAL'].values[0] / ceei_total

        #Append to final proportions dataframe
        to_proportions = pd.DataFrame(data={'REGION': [r], 'PROPORTION_RES': [proportion_res], 'PROPORTION_CSMI': [proportion_csmi]})
        proportions = pd.concat([proportions, to_proportions])
    
    proportions = proportions.set_index('REGION')

    #Comox Valley & Strathcona combine into Comox-Strathcona, Metro-Vancouver becomes GreaterVancouver (this is to line up with GADM naming convention)
    proportions.loc['Comox Valley'] += proportions.loc['Strathcona']
    proportions = proportions.drop('Strathcona')
    proportions = proportions.rename(index={'Comox Valley': 'Comox-Strathcona', 'Metro-Vancouver': 'GreaterVancouver'})

    #Quick fix for when proportions don't sum to 1 (currently it sums to ~0.993200, or 99.32%)
    proportions /= (proportions['PROPORTION_RES'].sum() + proportions['PROPORTION_CSMI'].sum())

    #Construct dataframes for hourly load. 27 subdivisions, each with residential and industrial loads
    regional_res = pd.DataFrame()
    regional_csmi = pd.DataFrame()

    #Apply the proportions to the hourly load data to get the disaggregated hourly load for each region in BC
    for region in proportions.index:
        p_res = proportions.loc[region]['PROPORTION_RES']
        p_csmi = proportions.loc[region]['PROPORTION_CSMI']

        #Spaces are removed in the column names, now all region names line up with the GADM names
        regional_res[region.replace(' ', '')] = hourly['LOAD'].apply(lambda x: x * p_res)
        regional_csmi[region.replace(' ', '')] = hourly['LOAD'].apply(lambda x: x * p_csmi)

    #Return the dataframes for disaggregated residential loads and disaggregated industrial loads
    return regional_res, regional_csmi


#Does some input verification before generating the regional hourly loads
'''

    USAGE:

    python dissagregate_load.py <ARG1> <ARG2> <ARG3> <ARG4>

    ARG1 = The CEEI spreadsheet file
    ARG2 = The folder which contains the hourly load data from BC Hydro. This folder should contain files named BalancingAuthorityLoad20XX.xls
    ARG3 = The year to use for disaggregation
    ARG4 = The folder to write the outputs to. This script outputs two CSV files, one for residential load, one for industrial load

    EXAMPLE:

    python disaggregate_load.py CEEI_2020.xlsx .\load\ 2015 .\hourly_load_region\
    
'''
def main():
    try:
        #Try reading the arguments passed in the terminal
        ceei_path = Path(sys.argv[1])
        hourly_path = Path(sys.argv[2] + '/BalancingAuthorityLoad' + sys.argv[3] + '.xls')
        year = int(sys.argv[3])
        output_path_res = Path(sys.argv[4] + '/hourly_res_' + sys.argv[3] + '.csv')
        output_path_csmi = Path(sys.argv[4] + '/hourly_csmi_' + sys.argv[3] + '.csv')

    except Exception as e:
        #Less than 4 arguments, return error code 1
        print('There are inputs missing')
        print(e)
        return 1
    
    else:
        #Correct number of arguments
        try:
            #Try loading in the data
            ceei = pd.read_excel(ceei_path, sheet_name='BC Hydro')
            ceei = ceei.loc[ceei.YEAR == min(2020, year)] #CEEI data only goes up to 2020, use 2020 if looking at hourly loads 2020 or later

            #Hourly load data needs some fixing
            hourly = fix_hourly_load(pd.read_excel(hourly_path), year)

        except Exception as e:
            #An input may be spelled incorrectly, return error code 2
            print('One or more inputs are in the wrong format')
            print(e)
            return 2
        else:
            #All is good, start hourly_res_{year}.csv and hourly_csmi_{year}.csv
            hourly_res, hourly_csmi = disaggregate(ceei, hourly)

            #Write to files to the output folder path
            hourly_res.to_csv(output_path_res)
            hourly_csmi.to_csv(output_path_csmi)
            
            #Return code 0 is for when everything runs without a problem
            return 0
        

if __name__ == '__main__':
    main()