import pandas as pd
"""
def load_and_process_data(file_path, 
                          technologies):
    df = pd.read_csv(file_path)
    year_col = 'YEAR'
    technology_col = 'TECHNOLOGY'
    if 'YEAR' not in df.columns or 'TECHNOLOGY' not in df.columns:
        year_col = 'y'
        technology_col = 't'
    yearly_summed_values = {tech: [] for tech in technologies}
    all_years = sorted(set(df[year_col]))
    grouped_data = {}
    for tech in technologies:
        filtered_df = df[df[technology_col].str.startswith(tech)]
        last_column = filtered_df.iloc[:, -1]
        grouped = last_column.groupby(filtered_df[year_col]).sum()
        grouped_data[tech] = grouped.reindex(all_years, fill_value=0)
        yearly_summed_values[tech] = [grouped.get(year, 0) for year in all_years]
    result_df = pd.DataFrame(grouped_data)

    return all_years, result_df

""" 

def load_and_process_data(file_path, 
                          technologies):
    df = pd.read_csv(file_path)
    year_col = 'YEAR'
    technology_col = 'TECHNOLOGY'
    if 'YEAR' not in df.columns or 'TECHNOLOGY' not in df.columns:
        year_col = 'y'
        technology_col = 't'
    yearly_summed_values = {tech: [] for tech in technologies}
    all_years = sorted(set(df[year_col]))
    grouped_data = {}
    for tech in technologies:
        filtered_df = df[df[technology_col].str.startswith(tech)]
        last_column = filtered_df.iloc[:, -1]
        grouped = last_column.groupby(filtered_df[year_col]).sum()
        grouped_data[tech] = grouped.reindex(all_years) #, fill_value=0
        yearly_summed_values[tech] = [grouped.get(year, 0) for year in all_years]
    result_df = pd.DataFrame(grouped_data)

    return all_years, result_df