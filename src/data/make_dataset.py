# -*- coding: utf-8 -*-
import click
import logging
import geopandas as gpd
from pathlib import Path
import pandas as pd
import os

#geodatabase for one of the data areas
ms4_gdb = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Tool_Preprocessing.gdb'
ms4_model_gdb = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'

print('Reading in municipal and regional boundaries...')

#from dotenv import find_dotenv, load_dotenv
munis_fp = "K:\\\DataServices\\Datasets\\Boundaries\\Spatial\\mapc_towns_poly.shp"
munis = gpd.read_file(munis_fp)

#mapc boundary 
mapc_boundary = gpd.read_file(ms4_gdb, layer='MAPCBoundary')

print('Reading in land parcel database data...')

#MAPC land parcel database csv
intermediate_path = "K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Data\\Intermediate"
mapc_lpd_fp = "K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Project_Files\\Parcel_Database\\MAPC_LandParcelDatabase.csv"
mapc_lpd = pd.read_csv(mapc_lpd_fp)


print('Reading in additional data layers...')
print('Wetlands...')
#wetlands (will need to pare down)
wetlands_fp = 'K:\\DataServices\\Datasets\\MassGIS\\Wetlands\\WETLANDSDEP_POLY.shp'
wetlands = gpd.read_file(wetlands_fp)
wetlands = wetlands.loc[wetlands['IT_VALDESC'] != 'OPEN WATER']
#add the clipping info at some point

print('Watersheds...')
#watersheds
major_basins_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Watersheds\\MAJBAS_POLY.shp'
major_basins = gpd.read_file(major_basins_fp)

subbasins_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Watersheds\\NRCSHUC10_POLY.shp'
subbasins = gpd.read_file(subbasins_fp)

print('Wellhead protection areas...')
#wellhead protection areas
interim_wpa_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Wellhead_Protection_Areas\\IWPA_POLY.shp'
zone1_wpa_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Wellhead_Protection_Areas\\ZONE1_POLY.shp'
zone2_wpa_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Wellhead_Protection_Areas\\ZONE2_POLY.shp'

interim_wpa = gpd.read_file(interim_wpa_fp)
zone1_wpa = gpd.read_file(zone1_wpa_fp)
zone2_wpa = gpd.read_file(zone2_wpa_fp)

#add this to make_dataset.py
interim_wpa['wpa_type'] = 'Interim Wellhead Protection Area'
zone1_wpa['wpa_type'] = 'Zone 1 Wellhead Protection Area'
zone2_wpa['wpa_type'] = 'Zone 2 Wellhead Protection Area'

combined_wpa = pd.concat([interim_wpa[['wpa_type', 'SITE_NAME', 'SUPPLIER', 'geometry']], 
                        zone1_wpa[['wpa_type', 'SITE_NAME', 'SUPPLIER','geometry']], 
                        zone2_wpa[['wpa_type', 'SUPPLIER','geometry']]]).reset_index()

print('Activity use limitation areas...')
#activity use limitation areas
aul_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\AUL\\AUL_PT.shp'
aul = gpd.read_file(aul_fp)

print('ParkServe data...')
#park serve priority areas data (this data is from 2022)
parkserve_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
parkserve_data = gpd.read_file(parkserve_fp, layer='TPL_parkserve_clip_mapc')

#do a tree need assessment
mapc_bgs_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
mapc_bgs = gpd.read_file(mapc_bgs_fp, layer='mapc_2020_blockgroups')

mapc_blocks_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
mapc_blocks = gpd.read_file(mapc_blocks_fp, layer='mapc_2020_blocks')

#environmental justice
ej_2020 = gpd.read_file(ms4_model_gdb, layer='ej_2020_bg_mapc')

#urban heat island index
heat_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Climate_Change\\MVP_MMC_Heat_MVP\\00 Task 2 Deliverables\\2.1 Attachments\\00 Uploaded to Sharepoint\\Shapefile_LSTIndex\\LSTindex.tif'
