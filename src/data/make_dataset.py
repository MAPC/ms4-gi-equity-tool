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

#lclu and soils layer within model gdb
lclu_layer = "lclu_simplify_all_mapc"
soils_layer = "soils_mapc_simplify"

#what field in lclu defines the phosphorus estimate?
pler_field = 'pler'

#what is the field that defines the hydrologic soil group?
hsg_field = 'HYDROLGRP'

#from lclu, what are the tree canopy cover names
tree_canopy_covernames = ['Deciduous Forest', 'Evergreen Forest', 'Palustrine Forested Wetland']

#right of way segments. these are derived from models > row_segmentation.py 
row_output_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\ROW_model_output.gdb'
row_segment_layer = 'ROW_segmentation_MAPC'
row_segment_layer_bos = 'ROW_segmentation_BOS' #boston has its own
row_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\RightOfWay_Segmentation.gdb'
eot_layer = 'EOTROADS_MAPC' #EOT roads layer

#drainage network geodatabase - layers were created manually in ArcGIS Pro
drainage_network_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\Drainage_network.gdb'

print('Reading in municipal and regional boundaries...')

#municipal boundaries for mapc region
munis_fp = "K:\\\DataServices\\Datasets\\Boundaries\\Spatial\\mapc_towns_poly.shp"
munis = gpd.read_file(munis_fp)

#mapc boundary 
mapc_boundary = gpd.read_file(ms4_gdb, layer='MAPCBoundary')

print('Reading in land parcel database data...')

#path to parcel data for section 3a. others are in ms4 gdb
section3a_parcels_path = "K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Data\\Intermediate"


print('Reading in additional data layers...')
print('Wetlands...')

#wetlands 
wetlands_fp = 'K:\\DataServices\\Datasets\\MassGIS\\Wetlands\\WETLANDSDEP_POLY.shp'
wetlands = gpd.read_file(wetlands_fp)
wetlands = wetlands.loc[wetlands['IT_VALDESC'] != 'OPEN WATER']
wetlands_field = 'IT_VALDESC'

print('Watersheds...')

#watershed
watersheds_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Watersheds\\NRCSHUC10_POLY.shp'
watersheds = gpd.read_file(watersheds_fp)
watersheds_field = 'HU_10_NAME'

print('Wellhead protection areas...')
#wellhead protection areas
interim_wpa_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Wellhead_Protection_Areas\\IWPA_POLY.shp'
zone1_wpa_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Wellhead_Protection_Areas\\ZONE1_POLY.shp'
zone2_wpa_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\Wellhead_Protection_Areas\\ZONE2_POLY.shp'

interim_wpa = gpd.read_file(interim_wpa_fp)
zone1_wpa = gpd.read_file(zone1_wpa_fp)
zone2_wpa = gpd.read_file(zone2_wpa_fp)

zone1_wpa_field = 'SUPPLIER'
zone2_wpa_field = 'SUPPLIER'
int_wpa_field = 'SUPPLIER'

print('Activity use limitation areas...')
#activity use limitation areas
aul_fp = 'K:\\DataServices\\Datasets\\Environment and Energy\\AUL\\AUL_PT.shp'
aul = gpd.read_file(aul_fp)
aul_field = 'NAME'

print('ParkServe data...')
#park serve priority areas data (this data is from 2022)
parkserve_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
parkserve_data = gpd.read_file(parkserve_fp, layer='TPL_parkserve_clip_mapc')
parkserve_field = 'ParkRank'

#do a tree need assessment
mapc_bgs_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
mapc_bgs = gpd.read_file(mapc_bgs_fp, layer='mapc_2020_blockgroups')

mapc_blocks_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
mapc_blocks = gpd.read_file(mapc_blocks_fp, layer='mapc_2020_blocks')

#environmental justice (2020 boundaries, updated in 2023)
ej_2020 = gpd.read_file(ms4_model_gdb, layer='ej_2020_bg_mapc')
ej_field = 'EJ_CRIT_DESC'

#urban heat island index
heat_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Climate_Change\\MVP_MMC_Heat_MVP\\00 Task 2 Deliverables\\2.1 Attachments\\00 Uploaded to Sharepoint\\Shapefile_LSTIndex\\LSTindex.tif'

#town centers
town_center_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Project_Files\\City_TownCenters\\city_towncenters.shp'
town_center = gpd.read_file(town_center_fp)

#community visibility layer, see R script src > features > ms4-comm-vis.R
community_vis_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Data\\Spatial\\Community_Visibility\\CommunityVisabilityLayer.shp'
community_vis = gpd.read_file(community_vis_fp)

#building structures
building_structures_gdb = 'K:\\DataServices\\Datasets\\MassGIS\\Facilities_Structures\\Building_Structures\\Output\\structures.gdb'
building_structures_layer = 'STRUCTURES_POLY'