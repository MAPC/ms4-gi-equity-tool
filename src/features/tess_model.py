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


from src.features.suitability_criteria import *
from src.data.make_dataset import *
from src.features.ms4_funcs import *


#see all columns in tables
pd.set_option('display.max_columns', None)



def tess_ms4_model(town_name, processed_path):
    '''
    
    description

    '''
    print('Starting on ' + town_name + '...')

    #create shapefile for muni
    muni_shp = munis.loc[munis['municipal'] == town_name]

    row_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\ROW_model_output.gdb'

    for layer_name in fiona.listlayers(row_gdb):
        if town_name in layer_name:
            town_row_layer = layer_name


    town_row = gpd.read_file(row_gdb, layer=town_row_layer)
    town_row = town_row.explode().reset_index(drop=True)

    ##FISHNET##

    fishnets_gdb = 'K:\\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Fishnets.gdb'

    for layer_name in fiona.listlayers(fishnets_gdb):
        if town_name in layer_name:
            town_tess = layer_name

    muni_tess = gpd.read_file(fishnets_gdb, layer=town_tess)
    #muni_tess = gpd.read_file(muni_tess_gdb, layer='muni_tessellation')
    #clip tesselation to the franklin boundary
    muni_tess = muni_tess.clip(muni_shp)
    muni_tess['sqm'] = muni_tess['geometry'].area

    ## PARCELS ##
    print('Prepping parcel data...')
    town_parcels_row = get_parcels_row(town_name, mapc_lpd)

    #create an imperviousness layer

    ## LAND COVER DATA PREP ##
    print('Prepping land cover data...')

    #first read in lclu with a muni shapefile mask
    lclu_muni = get_lclu_muni(muni_shp)

    #then get tree canopy only
    lclu_treecanopy = get_tree_canopy_lc(lclu_muni)

    #then calculate imperviousness 
    #select only for imperviousness land cover codes and fix gemoetry issues
    lclu_impervious = lclu_muni.loc[lclu_muni['covercode'] == 2]
    lclu_impervious.geometry = lclu_impervious.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

    #read in structures for the muni
    building_structures_gdb = 'K:\\DataServices\\Datasets\\MassGIS\\Facilities_Structures\\Building_Structures\\Output\\structures.gdb'
    building_structures = gpd.read_file(building_structures_gdb, layer='STRUCTURES_POLY', mask=muni_shp)

    #add a type field
    building_structures['type'] = 'rooftop'

    #erase structures from the imperviousness lclu layer
    #now we just have land cover (not rooftops) that are imperviousness
    imperv_cover = lclu_impervious.overlay(building_structures, how='difference')

    #add a type field
    imperv_cover['type'] = 'land cover'    


    print('Calculating parcel score for each criteria...')
    #does it overlap with public land?
    public = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='GRID_ID',
                                            parcel_data=muni_tess, 
                                            join_data=public_land, 
                                            field='USE_DESC',
                                            layer_name='pblc',
                                            overlap=0.05)
    town_row['ROW_ID'] = town_row['parloc_id']
    row_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                join_data=town_row, 
                                                field='ROW_ID',
                                                layer_name='row',
                                                overlap=0.05)

    #how much imperviousness is on site?
    imperv = calculate_imperviousness(muni_tess, 
                                    'GRID_ID', 
                                    imperv_cover, 
                                    building_structures)
    
    #what is the imperviousness pler on site?
    imperv_pler = get_K_load(parcel_data=muni_tess,
                            id_field = 'GRID_ID',
                            muni_name = town_name,
                            muni_shpf = muni_shp,
                            lclu_layer=lclu_impervious,
                            pler_field='pler')

    #does it overlap with Trust for Public Land's "Park Serve" priority areas from 2022
    parkserve = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='GRID_ID',
                                            parcel_data=muni_tess, 
                                            join_data=parkserve_data, 
                                            field='ParkRank',
                                            layer_name='parks',
                                            overlap=0.05)

    #does it overlap with wetlands?
    wetlands_ovlp = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                join_data=wetlands, 
                                                field='IT_VALDESC',
                                                layer_name='wtlnds',
                                                overlap=0.05)

    #does it overlap with wellhead protection areas?
    wpa_z1_ovlp_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                join_data=zone1_wpa, 
                                                field='SUPPLIER',
                                                layer_name='z1_wpa',
                                                overlap=0.05)

    wpa_z2_ovlp_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                join_data=zone2_wpa, 
                                                field='SUPPLIER',
                                                layer_name='z2_wpa',
                                                overlap=0.05)


    wpa_interim_ovlp_tess = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                join_data=interim_wpa, 
                                                field='SUPPLIER',
                                                layer_name='int_wpa',
                                                overlap=0.05)
    
    #does it have an activity use limitation area within it?
    aul_ovlp = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                points=True,
                                                join_data=aul, 
                                                field='NAME',
                                                layer_name='aul',
                                                overlap=0.05)

    #what watershed does it overlap with?
    watershed = calculate_suitability_criteria(how='if_overlap', 
                                                id_field='GRID_ID',
                                                parcel_data=muni_tess, 
                                                join_data=major_basins, 
                                                field='NAME',
                                                layer_name='wtshd',
                                                overlap=0.05)

    #is it in a block group with tree need?

    #get tree "need" layer (block group)
    muni_tree_need = tree_score(muni_shp, lclu_treecanopy)

    #calculate suitability
    tree_need = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='GRID_ID',
                                            parcel_data=muni_tess, 
                                            join_data=muni_tree_need, 
                                            field='rnk_tree',
                                            layer_name='tree',
                                            overlap=0.05)

    #is it in a census block with high heat vulnerability?

    #get heat "vulnerability" layer (block)
    muni_heat_vuln = heat_score(muni_shp, heat_fp)

    #calculate suitability
    heat_vln = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='GRID_ID',
                                            parcel_data=muni_tess, 
                                            join_data=muni_heat_vuln, 
                                            field='rnk_heat',
                                            layer_name='heat',
                                            overlap=0.05)
                                            
    #is it in an ej community?
    ej = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='GRID_ID',
                                            parcel_data=muni_tess, 
                                            join_data=ej_2020, 
                                            field='EJ_CRIT_DESC',
                                            layer_name='ej',
                                            overlap=0.05)
    
    #how close is it to the drainage network?
    #first identify the drainage network - let's just do lines for now.
    drainage_network_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\Drainage_network.gdb'
    for layer_name in fiona.listlayers(drainage_network_gdb):
        if town_name in layer_name:
            if 'lines' in layer_name:
                town_drainage_lines_layer = layer_name

    town_drainage_lines = gpd.read_file(drainage_network_gdb, layer=town_drainage_lines_layer)
    town_drainage_lines = town_drainage_lines.explode().reset_index(drop=True)

    drainage_lines = calculate_suitability_criteria(how='distance', 
                                                    id_field='GRID_ID',
                                                    parcel_data=muni_tess, 
                                                    join_data=town_drainage_lines,
                                                    layer_name='drn_ln')

    #how close is it to the nearest catch basin?
    for layer_name in fiona.listlayers(drainage_network_gdb):
        if town_name in layer_name:
            if 'pts' in layer_name:
                town_drainage_pts_layer = layer_name

    town_drainage_pts = gpd.read_file(drainage_network_gdb, layer=town_drainage_pts_layer)
    town_drainage_pts = town_drainage_pts.explode().reset_index(drop=True)

    drainage_pts = calculate_suitability_criteria(how='distance', 
                                                    id_field='GRID_ID',
                                                    parcel_data=muni_tess, 
                                                    join_data=town_drainage_pts,
                                                    layer_name='drn_pts')

    #how is soil quality?

    #get soil hydrology score
    muni_soils_hsg = soil_hsg_score(muni_shp)

    #then calculate
    soils = calculate_suitability_criteria(how='overlap_sjoin', 
                                            id_field='GRID_ID',
                                            parcel_data=muni_tess, 
                                            join_data=muni_soils_hsg, 
                                            field='hsg_scr',
                                            stats='mean',
                                            layer_name='soils',
                                            overlap=0.05)

    print('Building and exporting combined data table...')

    dfs = [row_tess,
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
        soils,
        drainage_lines,
        drainage_pts
        ]

    #merge all of the listed dataframes together based on their parcel id and geometry
    df_merged = reduce(lambda  left,right: pd.merge(left,
                                                    right,
                                                    on=['GRID_ID', 'geometry'], 
                                                    how='outer'), 
                                                    dfs)
    
    def conditions(row):
        if row['row'] == 1:
            val = 'ROW'
        elif row['pblc'] == 1:
            val = 'Public parcel'
        else:
            val = 'Private parcel'
        return val
        
    #add a "parcel type" field for pub, priv, and row
    df_merged['par_typ'] =  df_merged.apply(conditions, axis=1)

    df_merged.insert((len(df_merged.columns) - 1), 'geometry', df_merged.pop('geometry'))


    #export to shapefile
    path = os.path.join(processed_path, town_name) #make a subdirectory in intermediate folder w town name
    os.makedirs(path, exist_ok=True)
    df_merged.to_file(path + '\\' + town_name + '_MS4_GI_tesselation.shp')

    return df_merged
    