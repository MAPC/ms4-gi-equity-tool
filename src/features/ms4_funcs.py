
#build_features.py

import geopandas as gpd
import numpy as np
import pandas as pd
import os
import seaborn as sns
import pylusat #documentation: https://pylusat.readthedocs.io/en/latest/
import matplotlib.pyplot as plt
import geopandas as gpd
from scipy import stats
from scipy.stats import gamma
from pylusat import utils
import sklearn
import zipfile36 as zipfile
import fiona
import sys
sys.path.append("..")
from src import *
from functools import reduce
import rasterio
from rasterio.io import MemoryFile
import numpy as numpy
from rasterstats import zonal_stats
from rasterio import transform
from rasterio import features
from rasterio.enums import MergeAlg
from sklearn.preprocessing import MinMaxScaler
from shapely.validation import make_valid
from src.features.suitability_criteria import *

#from src.data.make_dataset import *

def get_parcels_row(town_name, muni_gdf):

    '''
    Creates base geometry for MS4 model. 

    INPUTS:
    town_name: string
    muni_gdf: geodataframe of municipal boundary

    OUTPUT:
    one geodataframe of parcels and right of way segments with
    unique identifiers ('parloc_id'), 'type' field that distinguishes parcels
    from ROW segments, and acreage. 
    
    For ROW segments, closest intersecting street name is added in 'address' field. 
    All fields retained for parcel data.

    '''

    ## RETRIEVE PARCEL DATA ##

    from src.data.make_dataset import section3a_parcels_path #this is where 3a parcel data lives
    from src.data.make_dataset import ms4_model_gdb 

    muni_id = muni_gdf.iloc[0]['muni_id'].astype(int) #get muni id to make row parcel id

    if town_name in os.listdir(section3a_parcels_path): #for towns in 3a list
        #prepare parcel layer wtih parcels and row segments
        town_parcels = get_landuse_data(town_name) #from src.features.build_features.py 
        #add muni and type fields
        town_parcels['muni'] = town_name 
        town_parcels['type'] = 'parcel'
    
    else: #for towns not in 3a list
        for layer_name in fiona.listlayers(ms4_model_gdb):
            if town_name in layer_name:
                town_parcels_layer = layer_name
        
        town_parcels = gpd.read_file(ms4_model_gdb, layer=town_parcels_layer)

        #rename field names to match fields in 3a parcel database
        town_parcels = town_parcels.rename(columns={'LOC_ID': 'parloc_id', 
                                                    'USE_DESC': 'UseDesc',
                                                    'TOWN_ID': 'muni_id',
                                                    'OWNER1': 'Owner',
                                                    'SITE_ADDR': 'Address'})
        town_parcels['muni'] = town_name
        town_parcels['type'] = 'parcel'
        town_parcels['acreage'] = town_parcels['geometry'].area / 4047
    

    ## RIGHT OF WAY ## 
    from src.data.make_dataset import row_output_gdb, row_segment_layer, row_segment_layer_bos, row_gdb, eot_layer
    
    #read in row segments, clip to municipal boundary
    if town_name == 'Boston': 
        print('getting row data from Boston')
        town_row = gpd.read_file(row_output_gdb, layer=row_segment_layer_bos, mask=muni_gdf)
    else:
        town_row = gpd.read_file(row_output_gdb, layer=row_segment_layer, mask=muni_gdf)
        town_row = town_row.drop(columns = 'parloc_id') #drop parlocid so you can recreate next
    town_row = town_row.clip(muni_gdf)

    #fix any geometry issues
    town_row = town_row.explode().reset_index(drop=True)
    town_row = town_row.loc[town_row['geometry'].geom_type == 'Polygon']

    #create ID field
    town_row.insert(0, 'parloc_id', range(10000, 10000 + len(town_row)))
    town_row['parloc_id'] = muni_id.astype(str) + '_ROW_' + town_row['parloc_id'].astype(str)
    town_row = town_row[['parloc_id', 'geometry']]

    #create type and muni fields that will be consistent with the parcels fields
    town_row['type'] = 'ROW segment'
    town_row['muni'] = town_name
    town_row['Owner'] = 'ROW segment'
    town_row['UseDesc'] = 'ROW segment'
    town_row['acreage'] = town_row['geometry'].area / 4047

    #join row segments to road EOT layer to get street name in address field
    eot_road_mapc = gpd.read_file(row_gdb, layer=eot_layer, mask=muni_gdf) #read in eot layer
    eot_road_mapc = eot_road_mapc.to_crs(town_parcels.crs)

    town_row = town_row.sjoin(eot_road_mapc, how="left")
    town_row['Address'] = town_row['STREETNAME']
    town_row = town_row[['parloc_id', 
                         'type', 
                         'muni', 
                         'Owner', 
                         'UseDesc',
                         'Address', 
                         'acreage', 
                         'geometry']].groupby(by='parloc_id').agg('first').reset_index()

    town_parcels = town_parcels.overlay(town_row, how='difference')

    #merge together ROW parcels with parcel data from 3a - final dataset for spatial operations
    town_parcels_row = pd.concat([town_parcels, town_row]).reset_index()
    town_parcels_row.geometry = town_parcels_row.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    
    return town_parcels_row

def get_lclu_muni(muni_gdf):
    '''
    for layer with 2016 lclu and pler estimates (created in models > pler_calc.py), 
    clips layer to municipal boundary and gives the total phosphorus load 
    '''
    from src.data.make_dataset import ms4_model_gdb, lclu_layer
    lclu = gpd.read_file(ms4_model_gdb, layer=lclu_layer, mask=muni_gdf)

    #fix any geometry issues
    lclu.geometry = lclu.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    #lclu = lclu.loc[lclu['geometry'].geom_type == 'Polygon']
    lclu['geometry'] = lclu['geometry'].buffer(0)
    lclu = lclu.clip(muni_gdf)
    #lclu.geometry = lclu.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    lclu = lclu.loc[lclu['geometry'].geom_type == 'Polygon']

    #get phosphorous load per land cover geometry area
    lclu['acres'] = lclu['geometry'].area / 4047
    lclu['P_load'] = lclu['pler'] * lclu['acres']
    return lclu

def get_tree_canopy_lc(lclu_muni):
    '''
    From municipal land cover dataset, extracts tree canopy land cover types. 
    Resulting geodataframe encompasses all area of land in 
    municipality that is covered by tree canopy.

    '''
    #tree canopy using land use
    
    from src.data.make_dataset import tree_canopy_covernames

    lclu_treecanopy = lclu_muni.loc[
        lclu_muni['covername'].isin(tree_canopy_covernames)]

    #fix geom inconsistnecies
    lclu_treecanopy = lclu_treecanopy.explode()
    lclu_treecanopy = lclu_treecanopy.loc[lclu_treecanopy['geometry'].geom_type == 'Polygon']
    lclu_treecanopy.geometry = lclu_treecanopy.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    return lclu_treecanopy


def tree_score(muni_boundary, muni_tree_canopy):
    '''
    
    For each block group in the municipality, determines the relative concentration of 
    tree canopy compared to all other block groups. Those in the bottom 40% of scores
    are retained as having tree "need". Parcels or fishnet cells within those block groups can
    then be prioritized higher.

    '''
    # PREPARE BASE GEOGRAPHY (BLOCK GROUPS) #

    #bring in block groups from dataset script
    from src.data.make_dataset import mapc_bgs

    # clip block groups to municipal boundary
    # need to eliminate "slivers" of block groups left from clip
    mapc_bgs['og_area'] = mapc_bgs['geometry'].area
    muni_bgs = mapc_bgs.clip(muni_boundary)
    muni_bgs['pct_bg'] = ((muni_bgs['geometry'].area) / (muni_bgs['og_area'])) * 100

    #only keep block groups where 5% or more of the bg remains
    muni_bgs = muni_bgs.loc[muni_bgs['pct_bg'] > 5]

    # TREE CANOPY COVERAGE #

    #use suitability calculation: for each block group, calculates a percentage of 
    #area that is covered by tree canopy


    muni_bgs_treecanopy = calculate_suitability_criteria(how='overlap_area', 
                                                        id_field='bg20_id',
                                                        parcel_data=muni_bgs, 
                                                        join_data=muni_tree_canopy, 
                                                        layer_name='tree')

    #do a percentile ranking and create categories based on the ranking
    muni_bgs_treecanopy['rnk_tree'] = muni_bgs_treecanopy['pct_tree'].rank(method='min', pct=True)

    tree_need_rule = [
                (muni_bgs_treecanopy['rnk_tree'] > 0.80),
                (muni_bgs_treecanopy['rnk_tree'] <= 0.80) & (muni_bgs_treecanopy['rnk_tree'] > 0.60),
                (muni_bgs_treecanopy['rnk_tree'] <= 0.60) & (muni_bgs_treecanopy['rnk_tree'] > 0.40),
                (muni_bgs_treecanopy['rnk_tree'] <= 0.40) & (muni_bgs_treecanopy['rnk_tree'] > 0.20),
                (muni_bgs_treecanopy['rnk_tree'] <= 0.20) & (muni_bgs_treecanopy['rnk_tree'] > 0)
        ]

    choices = ['Very high tree canopy', 'Moderately high tree canopy', 'Moderate tree canopy', 'Moderately low tree canopy', 'Very low tree canopy']

    muni_bgs_treecanopy['tree_need'] = np.select(tree_need_rule, choices, default=np.nan)

    #only keep bgs with the highest relative tree canopy need
    muni_bgs_tree_need = muni_bgs_treecanopy.loc[muni_bgs_treecanopy['rnk_tree'] <= 0.4]

    return muni_bgs_tree_need


def calculate_imperviousness(parcel_data, 
                             id_field, 
                             lclu_muni, 
                             muni_gdf):
    '''
    Inputs:
    - parcel data 
    - id_field for parcel data
    - lclu muni - lclu clipped to muni boundary
    - muni_gdf - boundary 
    - building structures in muni (pulled within function)

    Output - for each parcel, this function calculates:
    - Total area of impervious cover
    - Total area of impervious rooftop
    - Combined area of all impervious surface
    - Breakdown by percentage of impervious cover, rooftop, and pervious surface
    
    '''
    from src.data.make_dataset import building_structures_gdb, building_structures_layer

    #create impervious cover layer by exracting only impervious cover code (2)
    imprv_cover_layer = lclu_muni.loc[lclu_muni['covercode'] == 2]
    imprv_cover_layer.geometry = imprv_cover_layer.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

    #read in structures for the muni 
    building_structures = gpd.read_file(building_structures_gdb, 
                                        layer=building_structures_layer, 
                                        mask=muni_gdf)

    #add a type field
    building_structures['type'] = 'rooftop'

    #erase structures from the imperviousness lclu layer
    #now we just have land cover (not rooftops) that are imperviousness
    imperv_cover = imprv_cover_layer.overlay(building_structures, how='difference')

    #add a type field
    imperv_cover['type'] = 'land cover'  

    #calculate area of impervious cover
    imperv_surfaces = calculate_suitability_criteria(how='overlap_area', 
                                                id_field=id_field,
                                                parcel_data=parcel_data, 
                                                join_data=imprv_cover_layer, 
                                                field='type',
                                                layer_name='imp_cvr')

    #calculate area of rooftops
    imperv_rooftops = calculate_suitability_criteria(how='overlap_area', 
                                                    id_field=id_field,
                                                    parcel_data=parcel_data, 
                                                    join_data=building_structures, 
                                                    field='type',
                                                    layer_name='imp_rf')
    

    #join together based on parloc id
    imperv = imperv_surfaces.merge(imperv_rooftops[[id_field, 'sqm_imp_rf', 'pct_imp_rf']], 
                                on=id_field, 
                                how='outer').fillna(0)
        

    # use output from overlap_area and merge to create usable fields
    imperv['sqm_imprv'] = imperv['sqm_imp_cvr'] + imperv['sqm_imp_rf']
    imperv['pct_imprv'] = (imperv['sqm_imprv'] / imperv['geometry'].area) * 100
    imperv['sqm_prv'] = imperv['geometry'].area - imperv['sqm_imprv']
    imperv['pct_prv'] = 100 - imperv['pct_imprv']
    imperv['acr_imprv'] = imperv['sqm_imprv'] / 4047
    imperv['acr_imp_cvr'] = imperv['sqm_imp_cvr'] / 4047
    imperv['acr_imp_rf'] = imperv['sqm_imp_rf'] / 4047
    imperv['acr_prv'] = imperv['sqm_prv'] / 4047

    #move geometry to end
    imperv.insert((len(imperv.columns) - 1), 'geometry', imperv.pop('geometry'))

    return(imperv)

def heat_score (muni_boundary, heat_index_fp):

    '''
    For each census block  in the municipality, determines the relative heat index score compared to all other block groups. 
    Those in the top 40% of scores are retained as having heat "vulnerability". 
    Parcels within those block groups can then be prioritized higher.

    Inputs: Muni boundary (gdf), heat index raster (geotiff)

    '''
    from rasterstats import zonal_stats
    from src.data.make_dataset import mapc_blocks


    with rasterio.open(heat_index_fp) as raster:
        transform = raster.transform
        lst = raster.read(1).astype('float64')

    #read in census blocks, clip to muni and eliminate sliver blocks that remain
    mapc_blocks['og_area'] = mapc_blocks['geometry'].area
    muni_blocks = mapc_blocks.clip(muni_boundary)
    muni_blocks['pct_bg'] = ((muni_blocks['geometry'].area) / (muni_blocks['og_area'])) * 100

    #only keep block groups where 5% or more of the bg remains. Reset index for zonal stats
    muni_blocks = muni_blocks.loc[muni_blocks['pct_bg'] > 5].reset_index()

    #run zonal stats on heat index - what is the mean lst index score across census block?
    lst_stats = pd.DataFrame(zonal_stats(muni_blocks, 
                                        lst, 
                                        affine=transform, 
                                        stats='mean'))

    #join back to blocks, rename field
    muni_blocks_heat = muni_blocks[['geoid20', 'geometry']].join(lst_stats)
    muni_blocks_heat = muni_blocks_heat.rename(columns={'mean':'lst_mean'})

    #rank each block based on relative lst index score
    muni_blocks_heat['rnk_heat'] = muni_blocks_heat['lst_mean'].rank(method='min', pct=True)
    muni_blocks_heat.head()

    #create a categorical ranking for where each block lands relative to one another
    heat_rule = [
                    (muni_blocks_heat['rnk_heat'] > 0.80),
                    (muni_blocks_heat['rnk_heat'] <= 0.80) & (muni_blocks_heat['rnk_heat'] > 0.60),
                    (muni_blocks_heat['rnk_heat'] <= 0.60) & (muni_blocks_heat['rnk_heat'] > 0.40),
                    (muni_blocks_heat['rnk_heat'] <= 0.40) & (muni_blocks_heat['rnk_heat'] > 0.20),
                    (muni_blocks_heat['rnk_heat'] <= 0.20) & (muni_blocks_heat['rnk_heat'] > 0)
            ]

    choices = ['Very high heat inde', 'Moderately high heat index', 'Moderate heat index', 'Moderately low heat index', 'Very low heat index']

    muni_blocks_heat['heat_cmp'] = np.select(heat_rule, choices, default=np.nan)

    #only keep blocks with the highest relative heat index score
    muni_blocks_heat_vln = muni_blocks_heat.loc[muni_blocks_heat['rnk_heat'] >= 0.4]

    return(muni_blocks_heat_vln)


def soil_hsg_score(muni_shapefile):

    '''
    For each SSURGO geometry in the muni, calculates a soil hydrology
    score based on te soil hydrologic group. Soils in group A get a score of 3,
    soils in group B get a score of 2, soils in group C, D, A/D, B/D, and C/D
    get a score of 1. Soils without a hydrologic soil group get a score of 0.

    '''

    from src.data.make_dataset import ms4_gdb, soils_layer, hsg_field

    soils_hsg = gpd.read_file(ms4_gdb, layer=soils_layer, mask=muni_shapefile)

    soil_hsg_rule = [
        soils_hsg[hsg_field] == 'A',
        soils_hsg[hsg_field] == 'B',
        ((soils_hsg[hsg_field] == 'C') | 
        (soils_hsg[hsg_field] == 'D') | 
        (soils_hsg[hsg_field] == 'A/D') |
        (soils_hsg[hsg_field] == 'B/D') | 
        (soils_hsg[hsg_field] == 'C/D')),
        soils_hsg[hsg_field] == ' '
    ]

    choices = [3, 2, 1, 0]

    soils_hsg['hsg_scr'] = np.select(soil_hsg_rule, choices)

    return soils_hsg


def get_P_load (parcel_data, 
                id_field,
                muni_name,
                muni_gdf,
                lclu_layer):

    '''
    For each parcel, calculates phosphorous load using raster stats.

    Inputs:
    - Parcel data for muni + id field 
    - Name of muni (for naming clipped raster)
    
    Output: Input parcel data with id, geometry, and new fields:
    - P_per_acre - total Phosphorus load per acre of land on site
    - P_sum - total Phosphorus load on site

    Note that currently this is set to only calculate P load for impervious surfaces on site

    '''

    from rasterio.mask import mask
    from src.data.make_dataset import pler_field
    

    inRas = 'I:\\Elevation\\Lidar\\2023_MassGIS\\LiDAR_DEM_INT_16bit.tiff'
    outRas = 'I:\\Elevation\\Lidar\\2023_MassGIS\\' + muni_name + '_LiDAR_DEM_INT_16bit.tiff'


    #clip the dem to the muni, save to I drive folder
    with rasterio.open(inRas) as src:
        #update crs
        src.meta.update({
            'crs': muni_gdf.crs
            })
        out_image, out_transform= mask(src,muni_gdf.geometry,crop=True)
        out_meta=src.meta.copy() # copy the metadata of the source DEM

        
    out_meta.update({
        "driver":"Gtiff",
        "height":out_image.shape[1], # height starts with shape[1]
        "width":out_image.shape[2], # width starts with shape[2]
        "transform":out_transform
    })


    with rasterio.open(outRas,'w',**out_meta) as dst:
        dst.write(out_image)

    #need to transform pler into a raster otherwise this takes a really long time
    vector = lclu_layer[lclu_layer['covercode'] == 2]
    
    geom = [shapes for shapes in vector.geometry]
    geom_value = ((geom,value) for geom, value in zip(vector.geometry, vector[pler_field]))

    with rasterio.open(outRas) as raster:

        # Rasterize vector using the shape and transform of the raster
        rasterized = features.rasterize(geom_value,
                                        out_shape = raster.shape,
                                        transform = raster.transform,
                                        all_touched = True,
                                        fill = np.nan,   # background value
                                        merge_alg = MergeAlg.replace,
                                        dtype = np.float64)


        #run zonal statistics on pler raster
        pler_stats = pd.DataFrame(zonal_stats(parcel_data, 
                                                rasterized, 
                                                affine=raster.transform, 
                                                stats='sum'))

    parcel_pler_stats = parcel_data[[id_field, 'geometry']].join(pler_stats)
    parcel_pler_stats = parcel_pler_stats.rename(columns={'sum':'P_sum'})
    parcel_pler_stats['P_per_acre'] = parcel_pler_stats['P_sum'] / (parcel_pler_stats['geometry'].area / 4047)

    os.remove(outRas)

    return parcel_pler_stats



def comm_vis_layer(id_field, 
                   parcel_data, 
                   town_center_data,
                   comm_vis_data):

    '''
    For each parcel, calculates:
    - Whether or not the parcel is within a 'town center'
    - Calculates distance from nearest "visibile" community location 
    - If within threshold distance, provides information about what the nearest location is 

    Inputs:
    - Parcel data with ID field identified
    - Town center layer
    - Community visibility points layer (generated from src > features > ms4-comm-vis.R)
    
    Output - Parcel data with additional fields:
    - 'twncntr' - whether or not  (1, 0) site is within town center
    - 'comm' - whether or not (1,0) site is within 100 m from community asset
    - 'comm_name' and 'comm_type' - where 'comm' == 1, info about mearest community site
    - 'comm_viz' - whether or not (1,0) site is either within a town center OR within 100 m from community asset

    ''' 

    #distance to town centers

    towncenter = calculate_suitability_criteria(how='distance', 
                                                id_field=id_field,
                                                parcel_data=parcel_data, 
                                                join_data=town_center_data, 
                                                layer_name='twncntr')
    
    # create field for whether or not site is within town center (distance = 0) 
    towncenter_rule = [
                towncenter['dst_twncntr'] == 0,
                towncenter['dst_twncntr'] != 0
                ]
    choices = [1, 0]
    towncenter['twncntr'] = np.select(towncenter_rule, choices, default=np.nan)

    # define fields to retain in distance calculation
    comm_fields = ['NAME', 'TYPE']

    #calculate distance to nearest comm viz point, retain listed points
    comm_assets =  distance_with_fields (parcel_data = parcel_data, 
                                        id_field = id_field, 
                                        fields = comm_fields,
                                        distance_layer = comm_vis_data, 
                                        layer_name='comm')
    
    #determine threshold distance - currently set at 100m
    threshold_distance = 100

    # create field for whether site is within threshold distance from nearest comm. site
    comm_rule = [
                comm_assets['dst_comm'] <= threshold_distance,
                comm_assets['dst_comm'] > threshold_distance
                ]

    choices = [1, 0]
    comm_assets['comm'] = np.select(comm_rule, choices, default=np.nan)
    comm_assets = comm_assets.rename(columns = {"NAME": "comm_name", 
                                                "TYPE": "comm_type"}) #rename 

    # if within 100 m, retain name and type. otherwise leave null
    comm_assets['comm_name'] = np.where(comm_assets['comm']== 1, comm_assets['comm_name'], np.nan)
    comm_assets['comm_type'] = np.where(comm_assets['comm']== 1, comm_assets['comm_type'], np.nan)
    comm_vis = towncenter.merge(comm_assets, on=[id_field, 'geometry'], how='inner')

    # create a field for whether site is within town center OR within threshold distance from comm. site
    comm_vis_rule = [
        ((comm_vis['twncntr'] == 1) | (comm_vis['comm'] == 1)),
        ((comm_vis['twncntr'] != 1) & (comm_vis['comm'] != 1))
    ]

    choices = [1, 0]

    comm_vis['commvis'] = np.select(comm_vis_rule, choices, default=np.nan)
    
    return comm_vis

def drainage_data(town_name, 
                  id_field, 
                  parcel_data):
    '''
    Where available, reads in drainage data (points, lines) and calculates distance
    to drainage system for each site. If there is no drainage data available for the town,
    distance is set to -999.

    Inputs:
    - town name
    - id field
    - parcel data
    - drainage data (read in function, built manually)

    Output: input parcel data with additional fields:
    - 'dst_drn_ln': distance to nearest drainage pipe
    - 'dst_drn_pt': distance to nearest catch basin

    
    '''
    #first identify the drainage network 
    from src.data.make_dataset import drainage_network_gdb

    if any(town_name in s for s in fiona.listlayers(drainage_network_gdb)):
        #if town name is in the drainage network geodatabase, then run the following script to produce distance data

        #find layer in drainage geodatabase by looking for town name and 'lines'
        for layer_name in fiona.listlayers(drainage_network_gdb):
            if town_name in layer_name:
                if 'lines' in layer_name:
                    town_drainage_lines_layer = layer_name

        #read in  layer
        town_drainage_lines = gpd.read_file(drainage_network_gdb, layer=town_drainage_lines_layer)
        town_drainage_lines = town_drainage_lines.explode().reset_index(drop=True)

        #calculate distance from drainage lines 
        drainage_lines = calculate_suitability_criteria(how='distance', 
                                                        id_field=id_field,
                                                        parcel_data=parcel_data, 
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
                                                    id_field=id_field,
                                                    parcel_data=parcel_data, 
                                                    join_data=town_drainage_pts,
                                                    layer_name='drn_pts')
    else:
        drainage_lines = parcel_data[[id_field, 'geometry']].copy()
        drainage_lines['dst_drn_ln'] = -999
        drainage_lines['nrm_drn_ln'] = -999
        drainage_pts = parcel_data[[id_field, 'geometry']].copy()
        drainage_pts['dst_drn_pt'] = -999
        drainage_pts['nrm_drn_pt'] = -999

    drainage = drainage_lines.merge(drainage_pts, on=[id_field, 'geometry'], how='inner')
    drainage.insert((len(drainage.columns) - 1), 'geometry', drainage.pop('geometry'))

    return drainage