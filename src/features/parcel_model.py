'''
March to June 2023 
creator: rbowers

In support of the MS4 Equitable Green Infrastructure Site Selection Tool, this suitability model processes
statewide and local data to provide information about key characteristics related to green infrastructure. 

Geography: parcels and right of way segments

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
from datetime import datetime


from src.features.build_features import *
from src.data.make_dataset import *
from src.features.ms4_funcs import *


#see all columns in tables
pd.set_option('display.max_columns', None)

def parcel_ms4_model(town_name, processed_path):
    '''
    Parcel- and right-of-way based model that outputs information about various
    criteria related to equity other key characteristics for siting green infrastructure.
 
    Inputs:
    Town_name: string-based
    Processed_path: output path for where shapefile should be exported

    Output:
    Shapefile of parcels and right of way segments with fields for each tool criteria element.

    '''

    #starttime message
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print('Starting on ' + town_name + ' at ' + now + '...')

    #create geodataframe for muni boundary, used as input for several functions
    muni_shp = munis.loc[munis['municipal'] == town_name]

    ## PARCELS ## 
    town_parcels_row = get_parcels_row(town_name, muni_shp) #from src.features.ms4_funcs.py

    ## LAND COVER DATA PREP ##

    #first read in lclu with a muni shapefile mask
    lclu_muni = get_lclu_muni(muni_shp) #from src.features.ms4_funcs.py 

    #then get land cover for tree canopy only
    lclu_treecanopy = get_tree_canopy_lc(lclu_muni) #from src.features.ms4_funcs.py 

    # IMPERVIOUSNESS # 

    '''   
    #select only imperviousness land cover codes and fix gemoetry issues
    lclu_impervious = lclu_muni.loc[lclu_muni['covercode'] == 2]
    lclu_impervious.geometry = lclu_impervious.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

    #read in structures for the muni
    building_structures = gpd.read_file(building_structures_gdb, layer=building_structures_layer, mask=muni_shp)

    #add a type field
    building_structures['type'] = 'rooftop'

    #erase structures from the imperviousness lclu layer
    #now we just have land cover (not rooftops) that are imperviousness
    imperv_cover = lclu_impervious.overlay(building_structures, how='difference')

    #add a type field
    imperv_cover['type'] = 'land cover'    
    '''

    ## SUITABILITY MODEL ##

    #does it have a public land use? or a public owner signal?

    from src.data.public_uses import public_land_uses, owner_types

    town_parcels_row['pblc'] = town_parcels_row['UseDesc'].apply(lambda x: 1 if x in public_land_uses else 0).astype(int)
    town_parcels_row['pblc'] = town_parcels_row['Owner'].str.contains('|'.join(owner_types), na=False).astype(int)

    #how much imperviousness is on site?
    imperv = calculate_imperviousness(town_parcels_row, 
                                    'parloc_id', 
                                    lclu_muni, 
                                    muni_shp)
    
    #what is the imperviousness pler on site?
    imperv_pler = get_P_load(parcel_data=town_parcels_row,
                            id_field = 'parloc_id',
                            muni_name = town_name,
                            muni_gdf = muni_shp,
                            lclu_layer=lclu_muni,
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
                                                join_data=subbasins, 
                                                field='HU_10_NAME',
                                                layer_name='wtshd',
                                                overlap=0.05)

    #is it in a block group with tree need?

    #get tree "need" layer (block group)
    muni_tree_need = tree_score(muni_shp, lclu_treecanopy)

    #calculate suitability
    tree_need = calculate_suitability_criteria(how='if_overlap', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=muni_tree_need, 
                                            field='rnk_tree',
                                            layer_name='tree',
                                            overlap=0.05)

    #is it in a census block with high heat vulnerability?

    #get heat "vulnerability" layer (block)
    muni_heat_vuln = heat_score(muni_shp, heat_fp)

    #calculate suitability
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
    

    #is it visible in the community?
    comm_vis = comm_vis_layer(id_field = 'parloc_id', 
                                parcel_data = town_parcels_row, 
                                town_center_data = town_center,
                                comm_vis_data = community_vis)
    
    #how close is it to the drainage network?
    drainage = drainage_data(town_name=town_name, 
                            id_field='parloc_id', 
                            parcel_data=town_parcels_row)

    #how is soil quality?

    #get soil hydrology score
    muni_soils_hsg = soil_hsg_score(muni_shp)

    #then calculate the mean soil hydrology score across the parcel
    soils = calculate_suitability_criteria(how='overlap_sjoin', 
                                            id_field='parloc_id',
                                            parcel_data=town_parcels_row, 
                                            join_data=muni_soils_hsg, 
                                            field='hsg_scr',
                                            stats='mean',
                                            layer_name='soils',
                                            overlap=0.05)


    dfs = [town_parcels_row[['parloc_id', 
                             'Address', 
                             'Owner', 
                             'UseDesc', 
                             'muni', 
                             'type', 
                             'acreage', 
                             'pblc', 
                             'geometry']], 
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
        comm_vis,
        soils,
        drainage        
    ]

    #merge all of the listed dataframes together based on their parcel id and geometry
    df_merged = reduce(lambda  left,right: pd.merge(left,
                                                    right,
                                                    on=['parloc_id', 'geometry'], 
                                                    how='outer'), 
                                                    dfs)
    
    def conditions(row):
        if row['type'] == 'ROW segment':
            val = 'ROW segment'
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
    df_merged.to_file(path + '\\' + town_name + '_MS4_GI_parcels.shp')

    return df_merged



