
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
from src.features.build_features import *
#from src.data.make_dataset import *

def tree_score(muni_boundary, muni_tree_canopy):
    '''
    
    For each block group in the municipality, determines the relative concentration of 
    tree canopy compared to all other block groups. Those in the bottom 40% of scores
    are retained as having tree "need". Parcels or fishnet cells within those block groups can
    then be prioritized higher.

    '''
    #bring in block groups from dataset script
    from src.data.make_dataset import mapc_bgs

    mapc_bgs['og_area'] = mapc_bgs['geometry'].area
    muni_bgs = mapc_bgs.clip(muni_boundary)
    muni_bgs['pct_bg'] = ((muni_bgs['geometry'].area) / (muni_bgs['og_area'])) * 100

    #only keep block groups where 5% or more of the bg remains
    muni_bgs = muni_bgs.loc[muni_bgs['pct_bg'] > 5]

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


def calculate_imperviousness(parcel_data, id_field, imprv_cover_layer, imprv_structure_layer):
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
                                                    join_data=imprv_structure_layer, 
                                                    field='type',
                                                    layer_name='imp_rf')
    #join together based on parloc id
    imperv = imperv_surfaces.merge(imperv_rooftops[[id_field, 'sqm_imp_rf', 'pct_imp_rf']], 
                                on=id_field, 
                                how='inner')
    
    imperv['sqm_imprv'] = imperv['sqm_imp_cvr'] + imperv['sqm_imp_rf']
    imperv['pct_imprv'] = imperv['sqm_imprv'] / imperv['geometry'].area
    imperv['sqm_prv'] = imperv['geometry'].area - imperv['sqm_imprv']
    imperv['pct_prv'] = 1 - imperv['pct_imprv']
    imperv.insert((len(imperv.columns) - 1), 'geometry', imperv.pop('geometry'))
    return(imperv)

def heat_score (muni_boundary, heat_index_fp):
    '''
    
    For each block group in the municipality, determines the relative heat index score compared to all other block groups. 
    Those in the top 40% of scores
    are retained as having heat "vulnerability". Parcels or fishnet cells within those block groups can
    then be prioritized higher.

    '''
    from rasterstats import zonal_stats

    with rasterio.open(heat_index_fp) as raster:
        transform = raster.transform
        lst = raster.read(1).astype('float64')

    #read in census blocks, clip to muni and eliminate sliver blocks that remain
    mapc_blocks_fp = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
    mapc_blocks = gpd.read_file(mapc_blocks_fp, layer='mapc_2020_blocks')
    mapc_blocks['og_area'] = mapc_blocks['geometry'].area
    muni_blocks = mapc_blocks.clip(muni_boundary)
    muni_blocks['pct_bg'] = ((muni_blocks['geometry'].area) / (muni_blocks['og_area'])) * 100

    #only keep block groups where 5% or more of the bg remains. Reset index for zonal stats
    muni_blocks = muni_blocks.loc[muni_blocks['pct_bg'] > 5].reset_index()

    #run zonal stats on heat index
    lst_stats = pd.DataFrame(zonal_stats(muni_blocks, 
                                        lst, 
                                        affine=transform, 
                                        stats='mean'))

    #join back to blocks, rename field
    muni_blocks_heat = muni_blocks[['geoid20', 'geometry']].join(lst_stats)
    muni_blocks_heat = muni_blocks_heat.rename(columns={'mean':'lst_mean'})

    muni_blocks_heat['rnk_heat'] = muni_blocks_heat['lst_mean'].rank(method='min', pct=True)
    muni_blocks_heat.head()

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
    from src.data.make_dataset import ms4_gdb
    soils_hsg = gpd.read_file(ms4_gdb, layer="soils_mapc_simplify", mask=muni_shapefile)

    soil_hsg_rule = [
        soils_hsg['HYDROLGRP'] == 'A',
        soils_hsg['HYDROLGRP'] == 'B',
        ((soils_hsg['HYDROLGRP'] == 'C') | 
        (soils_hsg['HYDROLGRP'] == 'D') | 
        (soils_hsg['HYDROLGRP'] == 'A/D') |
        (soils_hsg['HYDROLGRP'] == 'B/D') | 
        (soils_hsg['HYDROLGRP'] == 'C/D')),
        soils_hsg['HYDROLGRP'] == ' '
    ]

    choices = [3, 2, 1, 0]

    soils_hsg['hsg_scr'] = np.select(soil_hsg_rule, choices)

    return soils_hsg


def get_K_load (parcel_data, 
                id_field,
                muni_name,
                muni_shpf,
                lclu_layer, 
                pler_field:str):

    '''
    describe here

    '''

    from rasterio.mask import mask
    

    inRas = 'I:\\Elevation\\Lidar\\2023_MassGIS\\LiDAR_DEM_INT_16bit.tiff'
    outRas = 'I:\\Elevation\\Lidar\\2023_MassGIS\\' + muni_name + '_LiDAR_DEM_INT_16bit.tiff'


    #clip the dem to the muni, save to I drive folder
    with rasterio.open(inRas) as src:
        #update crs
        src.meta.update({
            'crs': muni_shpf.crs
            })
        out_image, out_transform= mask(src,muni_shpf.geometry,crop=True)
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
    vector = lclu_layer
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
    parcel_pler_stats = parcel_pler_stats.rename(columns={'sum':'K_sum'})
    parcel_pler_stats['K_per_acre'] = parcel_pler_stats['K_sum'] / (parcel_pler_stats['geometry'].area / 4047)

    os.remove(outRas)

    return parcel_pler_stats
