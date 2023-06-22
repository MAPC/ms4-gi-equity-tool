
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
from src.features.build_features import *

#from src.data.make_dataset import *

def get_parcels_row(town_name, mapc_lpd, muni_shp):

    '''
    describe
    '''

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
    town_row = town_row[['poly_typ','geometry']]

    #create ID field
    muni_id = town_parcels.iloc[2]['muni_id'].astype(int)
    town_row.insert(0, 'parloc_id', range(10000, 10000 + len(town_row)))
    town_row['parloc_id'] = muni_id.astype(str) + '_ROW_' + town_row['parloc_id'].astype(str)

    #create type and muni fields that will be consistent with the parcels fields
    town_row['type'] = 'ROW segment'
    town_row['muni'] = town_name
    town_row['Owner'] = 'ROW segment'
    town_row['UseDesc'] = 'ROW segment'
    town_row['acreage'] = town_row['geometry'].area / 4047

    #join row segments to road EOT layer to get street name in address field
    row_gdb = 'K:\DataServices\Projects\Current_Projects\Environment\MS4\Project\RightOfWay_Segmentation.gdb'
    eot_road_mapc = gpd.read_file(row_gdb, layer='EOTROADS_MAPC', mask=muni_shp)
    eot_road_mapc = eot_road_mapc.to_crs(town_parcels.crs)

    town_row = town_row.sjoin(eot_road_mapc, how="left")
    town_row['Address'] = town_row['STREETNAME']
    town_row = town_row[['parloc_id', 'type', 'muni', 'Owner', 'Address', 'acreage', 'geometry']].groupby(by='parloc_id').agg('first').reset_index()


    #merge together ROW parcels with parcel data from 3a - final dataset for spatial operations
    town_parcels_row = pd.concat([town_parcels, town_row]).reset_index()
    town_parcels_row.geometry = town_parcels_row.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    
    return town_parcels_row

def get_lclu_muni(muni_gdf):
    '''
    describe
    '''
    ms4_model_gdb = 'K:\\DataServices\\Projects\\Current_Projects\\Environment\\MS4\\Project\\MS4_Model.gdb'
    lclu = gpd.read_file(ms4_model_gdb, layer="lclu_simplify_all_mapc", mask=muni_gdf)

    #fix any geometry issues
    lclu.geometry = lclu.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    lclu = lclu.clip(muni_gdf)
    lclu = lclu.loc[lclu['geometry'].geom_type == 'Polygon']

    #get phosphorous load per land cover geometry area
    lclu['acres'] = lclu['geometry'].area / 4047
    lclu['K_load'] = lclu['pler'] * lclu['acres']
    return lclu

def get_tree_canopy_lc(lclu_muni):
    '''
    describe
    '''
    #tree canopy using land use
    tree_canopy_covernames = ['Deciduous Forest', 'Evergreen Forest', 'Palustrine Forested Wetland']
    lclu_treecanopy = lclu_muni.loc[lclu_muni['covername'].isin(tree_canopy_covernames)]

    lclu_treecanopy.geometry = lclu_treecanopy.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)

    #fix geom inconsistnecies
    lclu_treecanopy = lclu_treecanopy.explode()
    lclu_treecanopy = lclu_treecanopy.loc[lclu_treecanopy['geometry'].geom_type == 'Polygon']
    return lclu_treecanopy


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
    '''
    describe
    '''
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
                                how='outer').fillna(0)
        

    
    imperv['sqm_imprv'] = imperv['sqm_imp_cvr'] + imperv['sqm_imp_rf']
    imperv['pct_imprv'] = (imperv['sqm_imprv'] / imperv['geometry'].area) * 100
    imperv['sqm_prv'] = imperv['geometry'].area - imperv['sqm_imprv']
    imperv['pct_prv'] = 100 - imperv['pct_imprv']
    imperv['acr_imprv'] = imperv['sqm_imprv'] / 4047
    imperv['acr_imp_cvr'] = imperv['sqm_imp_cvr'] / 4047
    imperv['acr_imp_rf'] = imperv['sqm_imp_rf'] / 4047
    imperv['acr_prv'] = imperv['sqm_prv'] / 4047

    imperv.insert((len(imperv.columns) - 1), 'geometry', imperv.pop('geometry'))

    return(imperv)

def heat_score (muni_boundary, heat_index_fp):
    '''

    For each block group in the municipality, determines the relative heat index score compared to all other block groups. 
    Those in the top 40% of scores
    are retained as having heat "vulnerability". Parcels or fishnet cells within those block groups can
    then be prioritized higher.

    Inputs: Muni boundary (gdf), heat index raster (geotiff)
    Outputs: 

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


public_land_uses = ['Vacant, Selectmen or City Council, Other City or Town (Municipal)',
    'Vacant, Selectmen or City Council (Municipal)',
    'Vacant, Conservation (Municipal or County)',
    'Improved, Selectmen or City Council (Municipal)',
    'Vacant, District (County)', 'United States Government',
    'Vacant, Education (Municipal or County)',
    'Improved, Education (Municipal or County)',
    'Improved, Other District (County)',
    'Vacant, Other District (County)',
    'Dept. of Conservation and Recreation (DCR) - Division of Water Supply Protection, Urban Parks (non-reimbursable)',
    'Mass. Highway Dept. (MHD) (non-reimbursable)',
    'Dept. of Conservation and Recreation (DCR) - Division of Urban Parks and Recreation (non-reimbursable)',
    'EXEMPT', 
    'Transportation Authority',
    'Improved, Municipal Public Safety', 'Improved, District (County)',
    'Dept. of Conservation and Recreation (DCR), Division of State Parks and Recreation',
    '(formerly Municipalities/Districts.  Removed June 2009.)',
    'Mass. Highway Dept. (MHD) (non-reimbursable), Gasoline Service Stations - providing engine repair or maintenance services, and fuel products',
    'Dept. of Fish and Game, Environmental Law Enforcement (DFG, formerly DFWELE) (non-reimbursable)',
    'Vacant, Selectmen or City Council (Municipal), Cemeteries (Charitable Org.)',
    'Recreation, Active Use (Charitable Org.)',
    'Improved, Selectmen or City Council, Other City or Town (Municipal)',
    'Non Profit Industrial', 
    'SEWER DEPT', 
    'UNKNOWN OWNER V',
    'municipal, Other Motor Vehicles Sales and Services',
    'TOWN-PROP  MDL-00',
    'Vacant, Conservation (Municipal or County), UNKNOWN OWNER V',
    'IMPUTED - Transportation Authority',
    'Vacant Land, UNKNOWN OWNER V',
    "Dept. of Corrections (DOC) - Division of Youth Services,Mass. Military,State Police,Sheriffs' Depts. (non-reimbursable)",
    'Utility Authority - Electric, Light, Sewer, Water',
    'Military Division - Campgrounds',
    'Comm. Of Mass. (Other, non-reimbursable)',
    'Developable Residential Land, (formerly Municipalities/Districts.  Removed June 2009.)',
    '(formerly Commonwealth of Massachusetts.  Removed June 2009.)',
    '(formerly Municipalities/Districts.  Removed June 2009.), Condo-Off',
    '(formerly Commonwealth of Massachusetts.  Removed June 2009.), Residential Condominium',
    'Single Family Residential, (formerly Municipalities/Districts.  Removed June 2009.)',
    'Condo-Off, (formerly Municipalities/Districts.  Removed June 2009.)',
    #added additional uses
    'Town Property Improved',
    'Bus Transportation Facilities and Related Properties',
    'IMPUTED - Housing Authority',
    'Housing Authority',
    'Vacant, Housing Authority',
    'Improved, Tax Title/Treasurer',
    'Vacant, Tax Title/Treasurer'
    ]

def comm_vis_layer(id_field, parcel_data, town_center_data,comm_vis_data):

    '''
    define this and add comments but this is the community visibility layer!
    inputs
    outputs

    ''' 

    #distance to town centers

    towncenter = calculate_suitability_criteria(how='distance', 
                                                id_field=id_field,
                                                parcel_data=parcel_data, 
                                                join_data=town_center_data, 
                                                layer_name='twncntr')

    towncenter_rule = [
                towncenter['dst_twncntr'] == 0,
                towncenter['dst_twncntr'] != 0
                ]

    choices = [1, 0]

    towncenter['twncntr'] = np.select(towncenter_rule, choices, default=np.nan)


    comm_fields = ['NAME', 'TYPE']

    comm_assets =  distance_with_fields (parcel_data = parcel_data, 
                                        id_field = id_field, 
                                        fields = comm_fields,
                                        distance_layer = comm_vis_data, 
                                        layer_name='comm')


    comm_rule = [
                comm_assets['dst_comm'] <= 100,
                comm_assets['dst_comm'] > 100
                ]

    choices = [1, 0]

    comm_assets['comm'] = np.select(comm_rule, choices, default=np.nan)
    comm_assets = comm_assets.rename(columns = {"NAME": "comm_name", "TYPE": "comm_type"})


    comm_vis = towncenter.merge(comm_assets, on=[id_field, 'geometry'], how='inner')

    comm_vis_rule = [
        ((comm_vis['twncntr'] == 1) | (comm_vis['comm'] == 1)),
        ((comm_vis['twncntr'] != 1) & (comm_vis['comm'] != 1))
    ]

    choices = [1, 0]

    comm_vis['commvis'] = np.select(comm_vis_rule, choices, default=np.nan)
    
    return comm_vis
