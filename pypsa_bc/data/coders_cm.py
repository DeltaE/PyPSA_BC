from pathlib import Path
import pandas as pd
import plotly.express as px
from shapely.geometry import Point
import geopandas as gpd
from linkingtool.coders import CODERSData
from linkingtool.AttributesParser import AttributesParser
from bc_combined_modelling.attributes_parser import AttributesParserExtended
from typing import Dict,Optional
import requests
from linkingtool.hdf5_handler import DataHandler
from dataclasses import dataclass, field

@dataclass
class CODERSData_extended(CODERSData,
                          AttributesParserExtended):
    
    def __init__(self, config_file_path,store_path:Optional[str] = None):
        
        self.config:Dict[str,dict] = self.load_config(config_file_path)
        
        if store_path is None:
            store_path= Path('data/downloaded_data.h5')
            print(f" >> Store not defined. Saving to default Local store :{store_path}")
            
        
        super().__post_init__()
        # self.data_pull_root = Path(self.data_pull['root'])
        self.store=DataHandler(store_path)
    

    
    def pull_all_data(self):
        
        for source, tables in self.data_pull.items():
            for table in tables:
                # Fetch the data using your existing methods
                data = self.fetch_data(table, source)
                print(f">> Retrieved data for {source}/{table} and saving to store")
                
                # Store the fetched data
                self.store.to_store(data, f'{source}/{table}',force_update=True)
         
    
    def fetch_data(self,
                   table,
                   source:str):
        """
        Args:
            table: table name
            source: cef, coders
        """
        if source=='coders':
            table_prefix=""
        else:
            table_prefix=f'/{source}'
        
        df = pd.DataFrame.from_dict(requests.get(self.url+f"{table_prefix}/{table}"+self.query).json())
        
        return df
    @staticmethod
    def convert_to_geodata(df):
        df['geometry'] = df.apply(lambda row: Point(row['longitude'], row['latitude']), axis=1)
        return gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326")

    def save_scenario_data(self, data: pd.DataFrame, table_name: str, scenario_col: str, folder: str):
        for scenario in data[scenario_col].unique():
            scenario_data = data[data[scenario_col] == scenario]
            file_path = self.data_pull_root / folder / f"{table_name}_{scenario}.csv"
            file_path.parent.mkdir(parents=True, exist_ok=True)
            scenario_data.to_csv(file_path, index=False)
            print(f"{table_name} scenario data saved to: {file_path}")

    def save_visualization(self, df: pd.DataFrame, x_col: str, y_cols: list, title: str, output_folder: str):
        fig = px.area(df, x=x_col, y=y_cols, title=title, labels={x_col: 'Year', 'value': 'Pj'})
        file_path = Path(output_folder) / f"{self.CRC}_visualization.html"
        file_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(file_path))
        print(f"Visualization saved to: {file_path}")


