'''
describe
'''


import sys
sys.path.append("..")
import os
os.environ['USE_PYGEOS'] = '0'
import pandas as pd
import geopandas as gpd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols 
from rasterstats import zonal_stats
import rasterio
import contextily as cx
from rasterio.plot import show
#import osmnx as ox
from affine import Affine
import rioxarray as rx
from multiprocessing import Pool, cpu_count
from shapely.ops import unary_union
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
from pathlib import Path
import fiona
from shapely.validation import make_valid


from src.features.build_features import *
from src.data.make_dataset import *
from src.features.ms4_funcs import *


#see all columns in tables
pd.set_option('display.max_columns', None)

def parcel_ms4_model(town_name):
    #create shapefile for muni
    print('Starting on ' + town_name + '...')
    muni_shp = munis.loc[munis['municipal'] == town_name]

    ## PARCELS ##
    print('Prepping parcel data...')

    row_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\ROW_model_output.gdb'

    for layer_name in fiona.listlayers(row_gdb):
        if town_name in layer_name:
            town_row_layer = layer_name

    town_row = gpd.read_file(row_gdb, layer=town_row_layer)
    town_row = town_row.explode().reset_index(drop=True)

    #prepare parcel layer wtih parcels and row segments
    town_parcels = get_landuse_data(town_name, mapc_lpd)
    town_parcels['muni'] = town_name
    town_parcels['type'] = 'parcel'


    town_parcels = town_parcels.overlay(town_row, how='difference')
    #set up row layer

    #don't need most field names
    town_row = town_row[['poly_typ', 'geometry']]

    #create ID field
    town_row.insert(0, 'parloc_id', range(10000, 10000 + len(town_row)))
    town_row['parloc_id'] = 'ROW_' + town_row['parloc_id'].astype(str)

    #create type and muni fields that will be consistent with the parcels fields
    town_row['type'] = 'ROW segment'
    town_row['muni'] = town_name

    #merge together ROW parcels with parcel data from 3a - final dataset for spatial operations
    town_parcels_row = pd.concat([town_parcels, town_row]).reset_index()
    town_parcels_row.geometry = town_parcels_row.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    #create an imperviousness layer

    ## LAND COVER DATA PREP ##
    print('Prepping land cover data...')
    #first read in lclu with a muni shapefile mask
    ms4_model_gdb = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
    lclu = gpd.read_file(ms4_model_gdb, layer="lclu_simplify_all_mapc", mask=muni_shp)

    #fix any geometry issues
    lclu.geometry = lclu.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    lclu = lclu.clip(muni_shp)
    lclu = lclu.loc[lclu['geometry'].geom_type == 'Polygon']

    #get phosphorous load per land cover geometry area
    lclu['acres'] = lclu['geometry'].area / 4047
    lclu['K_load'] = lclu['pler'] * lclu['acres']

    #select only for imperviousness land cover codes
    lclu_impervious = lclu.loc[lclu['covercode'] == 2]

    #fix any geometry issues
    lclu_impervious.geometry = lclu_impervious.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

    #tree canopy using land use
    tree_canopy_covernames = ['Deciduous Forest', 'Evergreen Forest', 'Palustrine Forested Wetland']
    lclu_treecanopy = lclu.loc[lclu['covername'].isin(tree_canopy_covernames)]

    lclu_treecanopy.geometry = lclu_treecanopy.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

    #fix geom inconsistnecies
    lclu_treecanopy = lclu_treecanopy.explode()
    lclu_treecanopy = lclu_treecanopy.loc[lclu_treecanopy['geometry'].geom_type == 'Polygon']

    #read in structures for the muni
    building_structures_gdb = 'K:\\DataServices\\Datasets\\MassGIS\\Facilities_Structures\\Building_Structures\\Output\\structures.gdb'
    building_structures = gpd.read_file(building_structures_gdb, layer='STRUCTURES_POLY', mask=muni_shp)

    #add a type field
    building_structures['type'] = 'rooftop'

    #erase structures from the imperviousness lclu layer
    imperv_cover = lclu_impervious.overlay(building_structures, how='difference')

    #now we just have land cover (not rooftops) that are imperviousness
    #add a type field
    imperv_cover['type'] = 'land cover'

    #join back together with the 'type' field being the distinguisher btwn land cover and rooftops
    imperviousness_by_type = pd.concat([building_structures[['type', 'geometry']],
                                        imperv_cover[['type', 'geometry']]
                                        ])

    #dissolve by type
    imperviousness_by_type = imperviousness_by_type.dissolve(by='type')

    #add a "layer" field for the model
    imperviousness_by_type['layer'] = 'impervious'

    print('Calculating relative tree canopy and heat vulnerability...')
    #get tree "need" layer (block group)
    muni_tree_need = tree_score(muni_shp, lclu_treecanopy)
    muni_tree_need.head()

    #get heat "vulnerability" layer (block)
    muni_heat_vuln = heat_score(muni_shp, heat_fp)

    #get soil hydrology score
    muni_soils_hsg = soil_hsg_score(muni_shp)

    print('Calculating parcel score for each criteria...')
    #does it overlap with public land?
    public = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=public_land, 
                                            field='USE_DESC',
                                            layer_name='pblc',
                                            overlap=0.05)
    

    #how much imperviousness is on site?
    imperv = calculate_imperviousness(town_parcels_row, 
                                    'parloc_id', 
                                    imperv_cover, 
                                    building_structures)
    #what is the imperviousness pler on site?
    imperv_pler = get_K_load(parcel_data=town_parcels_row,
                            id_field = 'parloc_id',
                            muni_name = town_name,
                            muni_shpf = muni_shp,
                            lclu_layer=lclu_impervious,
                            pler_field='pler')



    #does it overlap with Trust for Public Land's "Park Serve" priority areas from 2022
    parkserve = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=parkserve_data, 
                                            field='ParkRank',
                                            layer_name='parks',
                                            overlap=0.05)

    #does it overlap with wetlands?
    wetlands_ovlp = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='parloc_id',
                                                parcel_data=town_parcels_row, 
                                                join_data=wetlands, 
                                                field='IT_VALDESC',
                                                layer_name='wtlnds',
                                                overlap=0.05)

    #does it overlap with wellhead protection areas?
    wpa_z1_ovlp_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='parloc_id',
                                                parcel_data=town_parcels_row, 
                                                join_data=zone1_wpa, 
                                                field='SUPPLIER',
                                                layer_name='z1_wpa',
                                                overlap=0.05)

    wpa_z2_ovlp_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='parloc_id',
                                                parcel_data=town_parcels_row, 
                                                join_data=zone2_wpa, 
                                                field='SUPPLIER',
                                                layer_name='z2_wpa',
                                                overlap=0.05)


    wpa_interim_ovlp_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='parloc_id',
                                                parcel_data=town_parcels_row, 
                                                join_data=interim_wpa, 
                                                field='SUPPLIER',
                                                layer_name='int_wpa',
                                                overlap=0.05)
    #does it have an activity use limitation area within it?
    aul_ovlp = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='parloc_id',
                                                parcel_data=town_parcels_row, 
                                                points=True,
                                                join_data=aul, 
                                                field='NAME',
                                                layer_name='aul',
                                                overlap=0.05)

    #what watershed does it overlap with?
    watershed = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='parloc_id',
                                                parcel_data=town_parcels_row, 
                                                join_data=major_basins, 
                                                field='NAME',
                                                layer_name='wtshd',
                                                overlap=0.05)

    #is it in a block group with tree need?
    tree_need = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=muni_tree_need, 
                                            field='rnk_tree',
                                            layer_name='tree',
                                            overlap=0.05)

    #is it in a census block with high heat vulnerability?
    heat_vln = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=muni_heat_vuln, 
                                            field='rnk_heat',
                                            layer_name='heat',
                                            overlap=0.05)
                                            
    #is it in an ej community?
    ej = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=ej_2020, 
                                            field='EJ_CRIT_DESC',
                                            layer_name='ej',
                                            overlap=0.05)

    #how is soil quality?
    soils = calculate_suitability_criteria(how='overlap_sjoin', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=muni_soils_hsg, 
                                            field='hsg_scr',
                                            stats='mean',
                                            layer_name='soils',
                                            overlap=0.05)

    print('Building and exporting combined data table...')

    dfs = [town_parcels_row[['parloc_id', 'Address', 'Owner', 'UseDesc', 'muni', 'type', 'acreage', 'geometry']], 
        public,
        imperv,
        imperv_pler,
        parkserve,
        wetlands_ovlp,
        wpa_z1_ovlp_tess,
        wpa_z2_ovlp_tess,
        wpa_interim_ovlp_tess,
        aul_ovlp,
        watershed,
        tree_need,
        heat_vln,
        ej,
        soils
        ]

    #merge all of the listed dataframes together based on their parcel id and geometry
    df_merged = reduce(lambda  left,right: pd.merge(left,
                                                    right,
                                                    on=['parloc_id', 'geometry'], 
                                                    how='outer'), 
                                                    dfs)

    df_merged.insert((len(df_merged.columns) - 1), 'geometry', df_merged.pop('geometry'))

    #run suitability model for each town in list
    processed_path = "K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Data\\Spatial\\Output"

    #add back merged shapefile (add this back later)

    #export to shapefile
    path = os.path.join(processed_path, town_name) #make a subdirectory in intermediate folder w town name
    os.makedirs(path, exist_ok=True)
    df_merged.to_file(path + '\\' + town_name + '_MS4_GI_parcels.shp')

    return df_merged
    