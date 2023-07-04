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
from rasterio import transform
from rasterio import features
from rasterio.enums import MergeAlg

from sklearn.preprocessing import MinMaxScaler
#from src.data.make_dataset import *

def get_landuse_data(muni):
    '''
    input = muni name
    output = - picks out the right shapefile from the state's municipal land use database (for 3a);
             - makes a subdirectory in seciton 3a parcels folder w town name and exports land use shapefile to it
             - reads that shapefile in as a geodataframe

    '''


    from src.data.make_dataset import section3a_parcels_path

    path = os.path.join(section3a_parcels_path, muni)

    #set land  use variable
    town_lu_path = ""
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if filename.endswith('.shp'):
                town_lu_path = os.path.join(dirpath, filename)

    muni_state_parcels = gpd.read_file(town_lu_path).rename(columns={'LOC_ID': 'parloc_id'})
    
    return(muni_state_parcels)  
    

def rasterize_geom(vector, value_field, raster):


    #create tuples of geometry, value pairs, where value is the attribute value you want to burn
    #vector = lclu
    #value_field = 'COVERCODE'
    #raster = clipped

    geom_value = ((geom,value) for geom, value in zip(vector.geometry, vector[value_field]))
    
    rasterized_geom = features.rasterize(geom_value,
                                        out_shape = raster.shape,
                                        transform = raster.transform,
                                        all_touched = True,
                                        fill = -5,   # background value
                                        merge_alg = MergeAlg.replace,
                                        dtype = np.float32)
    
    return rasterized_geom

def overlap_sjoin_hold (target_layer, 
                   overlap_layer, 
                   field:str, 
                   stats=None):


    '''
    assigns the parcel with the (max, min, average, exact) value from an overlapping geography that may have multiple values overlapping the parcel.
    Best suited for an overlapping geography that may have multiple values within the parcel.

    inputs: 
    - target_layer = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - overlap_layer = the layer you are interested in getting the overlap (ie, transit buffer, flood zone, etc)
    - field (str) = which field you want to do statistics for overlap 
    - layer_name = name for the indicator
    - stats (str) = {'first', 'last', 'sum', 'mean', 'median', 'max', 'min', 'std', 'var', 'count', 'size'}. Default is none. 

    '''
    #from pylusat import geotools
    from src.features.sjoin_update import spatial_join

    if stats is None:
        stats = 'mean'
    else:
        stats = stats


    #valid = {'first', 'last', 'sum', 'mean', 'median', 'max', 'min', 'std', 'var', 'count', 'size'}
    #if stats not in valid:
    #    raise ValueError("stats must be one of %r." % valid)
        

    #reproject all to mass mainland
    mass_mainland_crs = "EPSG:26986"
    target_layer = target_layer.to_crs(mass_mainland_crs)
    overlap_layer = overlap_layer.to_crs(mass_mainland_crs)

    #only keep parts of overlap layer that intersects with parcels. using pylusat methodology but udpated
    overlap_geo_overlay = spatial_join(target_layer, 
                                        overlap_layer, 
                                        op='intersects', 
                                        cols_agg={field: [stats]},   
                                        join_type='one to one', 
                                        keep_all=True)
   
    

    return overlap_geo_overlay



def scaling(df, col:str):
    
    '''
    removes outliers then rescales column to a value from 0-1

    input = data frame with column name
    output = normalized value btwn 0 and 1 

    we can play around with methods layer. for now, it's min max scaling 
    https://towardsdatascience.com/data-normalization-with-pandas-and-scikit-learn-7c1cc6ed6475
    '''
    #cap outliers at Q1 - 1.5*IQR and Q3 + 1.5*IQR
    Q1=df[col].quantile(0.25) 
    Q3=df[col].quantile(0.75)
    IQR=Q3-Q1

    low_limit=Q1-1.5*IQR
    high_limit=Q3+1.5*IQR

    def trim_outliers (row):
        if row[col] < low_limit:
            return low_limit
        elif row[col] > high_limit:
            return high_limit
        else:
            return row[col]
        
    df_norm = df.copy()

    df_norm[col] = df_norm.apply(lambda row: trim_outliers(row), axis=1)
  
    #df_norm=df[col][~((df[col]<(Q1-1.5*IQR)) | (df[col]>(Q3+1.5*IQR)))]

    
    # apply min-max scaling to capped values
    df_norm = (df_norm[col] - df_norm[col].min()) / (df_norm[col].max() - df_norm[col].min())
        
    return df_norm



def calculate_overlap (parcel_data,
                       id_field, 
                       overlap_layer, 
                       layer_name:str, 
                       inverse=False):

    '''
    calculates the total area of overlap and the percentage of overlap for any polygon
    layer with a parcel layer.

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - overlap_layer = the layer you are interested in getting the overlap (ie, transit buffer, flood zone, etc)
    - layer_name (str) = for field names, a short name to indicate details about the overlay layer

    output:
    data frame with parcel id, area of parcel overlapped by overlay layer, percent of parcel overlapped by overlay layer, and 
    percentile ranking for area of parcel overlapped by overlay layer. 

    example = calculate percent of parcel that is within the 1/2 mile transit radius buffer

    '''
   
    #reproject all to mass mainland
    mass_mainland_crs = "EPSG:26986"
    parcel_data = parcel_data.to_crs(mass_mainland_crs)
    overlap_layer = overlap_layer.to_crs(mass_mainland_crs)

    #only keep parts of overlap layer that intersects with parcels
    overlap_parcels_overlay = parcel_data.overlay(overlap_layer, how='intersection')
    
    #get area of overlap for each parcel
    overlap_parcels_overlay['sqm_' + layer_name] = overlap_parcels_overlay['geometry'].area  
    overlap_parcels_overlay = overlap_parcels_overlay.groupby(by=id_field).agg({('sqm_' + layer_name):'sum'}).reset_index()
    

    #join back to parcels data, remove additional rows with groupby
    parcels_with_overlap = parcel_data.merge(overlap_parcels_overlay[[id_field, ('sqm_' + layer_name)]], on=id_field, how='left')

    #get percent of overlap for each aprcel
    parcels_with_overlap['pct_' + layer_name] = (parcels_with_overlap['sqm_' + layer_name] / (parcels_with_overlap['geometry'].area)) * 100

    #define final table
    parcels_with_overlap = parcels_with_overlap[[id_field, ('sqm_' + layer_name), ('pct_' + layer_name), 'geometry']]

    return parcels_with_overlap



def overlap_sjoin (parcel_data, 
                   id_field,
                   overlap_layer, 
                   field:str, 
                   layer_name:str, 
                   stats=None):

    '''
    assigns the parcel with the (max, min, average, exact) value from an overlapping geography that may have multiple values overlapping the parcel.
    Best suited for an overlapping geography that may have multiple values within the parcel.

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - overlap_layer = the layer you are interested in getting the overlap (ie, transit buffer, flood zone, etc)
    - field (str) = which field you want to do statistics for overlap 
    - layer_name = name for the indicator
    - stats (str) = {'first', 'last', 'sum', 'mean', 'median', 'max', 'min', 'std', 'var', 'count', 'size'}. Default is none. 

    '''
    from pylusat import geotools
    from src.features.sjoin_update import spatial_join

    if stats is None:
        stats = 'mean'
    else:
        stats = stats


    #valid = {'first', 'last', 'sum', 'mean', 'median', 'max', 'min', 'std', 'var', 'count', 'size'}
    #if stats not in valid:
    #raise ValueError("stats must be one of %r." % valid)
        

    #reproject all to mass mainland
    mass_mainland_crs = "EPSG:26986"
    parcel_data = parcel_data.to_crs(mass_mainland_crs)
    overlap_layer = overlap_layer.to_crs(mass_mainland_crs)

    #only keep parts of overlap layer that intersects with parcels. using pylusat methodology but udpated
    overlap_parcels_overlay = spatial_join(parcel_data, 
                                            overlap_layer, 
                                            op='intersects', 
                                            cols_agg={field: [stats]},   
                                            join_type='one to one', 
                                            keep_all=True)
    
    parcels_with_overlap = parcel_data.merge(overlap_parcels_overlay[[id_field, (field + '_' + stats)]], on=id_field, how='left').fillna(0)

    #get normalized value for each parcel 
    #parcels_with_overlap['nrm_' + layer_name] = scaling(parcels_with_overlap, (field + '_' + stats))

    #define final table
    parcels_with_overlap = parcels_with_overlap[[id_field, (field + '_' + stats), 'geometry']].fillna(0)

    return parcels_with_overlap
    

def distance (parcel_data, 
              id_field, 
              distance_layer, 
              layer_name:str):

    '''
    caclulates the distance from parcel boundaries to the nearest boundary of distance_layer

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - distance_layer = the layer you are interested in getting the distance to

    output:
    - distance field (miles) for distance in miles from parcel boundary to nearest distance_layer boundary
    - (inverse) normalized field. values closer to 1 = parcels closest to nearest distance_layer feature. values closer to 0 = parcels furthest
    from nearest distance_layer feature

    same as "near" in ArcGIS Pro https://pro.arcgis.com/en/pro-app/latest/tool-reference/analysis/near.htm

    '''
    #reproject all to mass mainland
    mass_mainland_crs = "EPSG:26986"
    parcel_data = parcel_data.to_crs(mass_mainland_crs)
    distance_layer = distance_layer.to_crs(mass_mainland_crs)

    #make a copy
    parcel_distance = parcel_data.copy()

    #use geopandas sjoin to get distance to nearest distance_layer
    
    parcel_distance = gpd.sjoin_nearest(parcel_data, distance_layer, how='left', distance_col = ('dst_' + layer_name))
    parcel_distance = parcel_distance.groupby(parcel_distance.index).agg('first')

    parcel_distance['dst_' + layer_name] = parcel_distance['dst_' + layer_name] 

    #normalize distance. want to prioritize closer distances so subtract from 1 
    parcel_distance['nrm_' + layer_name] = 1 - scaling(parcel_distance, ('dst_'+ layer_name))

    parcel_distance = parcel_data.merge(parcel_distance[[id_field, ('dst_' + layer_name), ('nrm_' + layer_name)]], on=id_field, how='left').fillna(0)

    #define final table
    parcel_distance = parcel_distance[[id_field, ('dst_'+ layer_name), ('nrm_' + layer_name), 'geometry']]

    return parcel_distance


def distance_with_fields (parcel_data, 
                            id_field, 
                            fields:list,
                            distance_layer, 
                            layer_name:str):

    '''
    caclulates the distance from parcel boundaries to the nearest boundary of distance_layer

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - distance_layer = the layer you are interested in getting the distance to
    - fields(list) = list of fields you want to retain from the nearest geometry in the distance_layer 

    output:
    - distance field (meters)) for distance in meters from parcel boundary to nearest distance_layer boundary
    - desired information about nearest distance_layer geometry

    from nearest distance_layer feature

    same as "near" in ArcGIS Pro https://pro.arcgis.com/en/pro-app/latest/tool-reference/analysis/near.htm

    '''

    #reproject all to mass mainland
    mass_mainland_crs = "EPSG:26986"
    parcel_data = parcel_data.to_crs(mass_mainland_crs)
    distance_layer = distance_layer.to_crs(mass_mainland_crs)

    #make a copy
    parcel_distance = parcel_data.copy()

    #use geopandas sjoin to get distance to nearest distance_layer
    
    parcel_distance = gpd.sjoin_nearest(parcel_data, distance_layer, how='left', distance_col = ('dst_' + layer_name))
    parcel_distance = parcel_distance.groupby(parcel_distance.index).agg('first')

    parcel_distance['dst_' + layer_name] = parcel_distance['dst_' + layer_name] 

    fields_list = [id_field, ('dst_' + layer_name)]
    
    for field in fields:
        fields_list.append(field)

    parcel_distance = parcel_data.merge(parcel_distance[fields_list], on=id_field, how='left').fillna(0)

    fields_list.append('geometry')

    #define final table
    parcel_distance = parcel_distance[fields_list]

    return parcel_distance

def parcel_ptile_table_outliers(parcel_data, 
                                id_field,
                                field_name:str, 
                                layer_name:str, 
                                inverse=False, 
                                nan_value=False):
    '''
    caclulates the percentile rank for a field within the initial parcel database. best for fields with outliers.

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - field_name = the field name to be normalized
    - layer_name = select a name for the normalized field. can be more read-able than 
    - inverse(binary) = if true, normalizes the inverse of the field so that normalized values closer to 0 are for higher values in the field. 

    output:
    - initial field
    - percentile rank
    '''

    #make a copy
    if nan_value:
        parcel_copy = parcel_data.replace(nan_value, np.NaN).copy()
        parcel_copy = parcel_copy.loc[parcel_copy[field_name] != np.NaN] #do percentile ranking without null values
    else:
        parcel_copy = parcel_data.copy()

    if inverse:
        parcel_copy['nrm_' + layer_name] = 1 - scaling(parcel_copy, field_name)
    else:
        parcel_copy['nrm_' + layer_name] = scaling(parcel_copy, field_name)

    '''if inverse:
        parcel_copy['ptl_' + layer_name] = 1 - parcel_copy[field_name].rank(method='min', pct=True)
    else:
        parcel_copy['ptl_' + layer_name] = parcel_copy[field_name].rank(method='min', pct=True)'''

    #merge to original table
    parcels_with_ptile = parcel_data.merge(parcel_copy[[id_field, ('nrm_' + layer_name)]], on=id_field, how='left')

    #define final table
    parcels_with_ptile = parcels_with_ptile[[id_field, field_name, ('nrm_' + layer_name), 'geometry']]

    return parcels_with_ptile

def parcel_norm_table(parcel_data, 
                      id_field,
                      field_name:str, 
                      layer_name:str, 
                      inverse=False, 
                      nan_value=False):
    '''
    caclulates the normalized value for a field within the initial parcel database

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - field_name = the field name to be normalized
    - layer_name = select a name for the normalized field. can be more read-able than 
    - inverse(binary) = if true, normalizes the inverse of the field so that normalized values closer to 0 are for higher values in the field. 

    output:
    - initial field
    - normalized field (using robust scaling)
    '''

    #make a copy
    if nan_value:
        parcel_copy = parcel_data.replace(nan_value, np.NaN).copy()
        parcel_copy = parcel_copy.loc[parcel_copy[field_name] != np.NaN] #do normalization without null values
    else:
        parcel_copy = parcel_data.copy()

    if inverse:
        parcel_copy['nrm_' + layer_name] = 1 - scaling(parcel_copy, field_name)
    else:
        parcel_copy['nrm_' + layer_name] = scaling(parcel_copy, field_name)

   
    #merge to original table
    parcels_with_ptile = parcel_data.merge(parcel_copy[[id_field, ('nrm_' + layer_name)]], on=id_field, how='left').fillna(0)


    #define final table
    parcels_with_ptile = parcels_with_ptile[[id_field, field_name, ('nrm_' + layer_name), 'geometry']]

    return parcels_with_ptile

def if_overlap (parcel_data, 
                id_field,
                overlap_layer, 
                field:str, 
                layer_name:str, 
                overlap=None, 
                points=False):

    '''
    assigns the parcel with a value of 1 if it overlaps with a layer of choice

    inputs: 
    - parcel_data = parcel database that you are joining overlap layer with (can be pre-filtered for muni of interest for speed)
    - overlap_layer = the layer you are interested in getting the overlap (ie, transit buffer, flood zone, historic points, etc)
    - points = define if the overlap field is a points layer or not (default = False)
    - field (str) = which field you want to get info from
    - layer_name (str) = for field names, a short name to indicate details about the overlay layer
    - overlap = define from 0 - 1, how much overlap is defined as "overlap" 

    '''
   
    #reproject all to mass mainland
    mass_mainland_crs = "EPSG:26986"
    parcel_data = parcel_data.to_crs(mass_mainland_crs)
    overlap_layer = overlap_layer.to_crs(mass_mainland_crs)

    valid = {
        True,
        False
        }
    
    if points not in valid:
        raise ValueError("points must be one of %r." % valid)
    
    if not points:
        if overlap is None:
            overlap = 0.5
        else:
            overlap = overlap

        #only keep parts of overlap layer that intersects with parcels
        overlap_parcels_overlay = parcel_data.overlay(overlap_layer, how='intersection')
        
        #calculate area of parcel overlap
        overlap_parcels_overlay['sqm_' + layer_name] = overlap_parcels_overlay['geometry'].area  
        overlap_parcels_overlay = overlap_parcels_overlay.groupby(by=id_field).agg({('sqm_' + layer_name):'sum', field:'first'}).reset_index()
    
        #join back to parcels data
        parcels_with_overlap = parcel_data.merge(overlap_parcels_overlay[[id_field, ('sqm_' + layer_name), field]], on=id_field, how='left').fillna(0)

        #get percent of overlap for each parcel
        parcels_with_overlap['pct_' + layer_name] = (parcels_with_overlap['sqm_' + layer_name] / (parcels_with_overlap['geometry'].area)) 

        #assign a value of 1 for parcels that overlap more than 50% with overlay layer
        parcels_with_overlap[layer_name] = np.where(parcels_with_overlap['pct_' + layer_name] > overlap, 1, 0).astype(int)   

        #define final table
        parcels_with_overlap = parcels_with_overlap[[id_field, field, layer_name, 'geometry']]

        return parcels_with_overlap
    
    else:
        point_in_poly = gpd.sjoin(parcel_data, overlap_layer, predicate='contains', how='inner')
        point_in_poly = point_in_poly.groupby(by=id_field).agg('first').reset_index()
    
        #join back to parcels data
        parcels_with_overlap = parcel_data.merge(point_in_poly[[id_field, field]], on=id_field, how='left')
        

        #assign a value of 1 for parcels that overlap with pt data
        parcels_with_overlap[layer_name] = np.where(parcels_with_overlap[field].isna(), 0, 1)   

        #define final table
        parcels_with_overlap = parcels_with_overlap[[id_field, layer_name, field, 'geometry']]

        return parcels_with_overlap



def calculate_suitability_criteria(
                                how,
                                parcel_data,  
                                layer_name:str, 
                                id_field,
                                join_data=None,
                                field=None, 
                                points=False,
                                stats=None,
                                inverse=None,
                                nan_value=None,
                                overlap=None
                                ):
    '''
    combines all of the different suitability criteria calculation
    types into one single function. 

    how:
    - 'overlap_area'
    - 'overlap_sjoin'
    - 'distance'
    - 'parcel_ptile_rank'
    - 'parcel_norm'
    - 'if_overlap'
    
    '''

    valid = {
        'overlap_area', 
        'overlap_sjoin', 
        'distance', 
        'parcel_ptile_rank', 
        'parcel_norm', 
        'if_overlap'
        }
    
    if how not in valid:
        raise ValueError("how must be one of %r." % valid)
    
    if how == 'overlap_area':
        return calculate_overlap(parcel_data, id_field, join_data, layer_name, inverse)
    
    elif how == 'overlap_sjoin':
        return overlap_sjoin(parcel_data=parcel_data, id_field=id_field,  overlap_layer=join_data, field=field, stats=stats, layer_name=layer_name)
    
    elif how == 'distance':
        return distance(parcel_data, id_field, join_data, layer_name)
    
    elif how == 'parcel_ptile_rank':
        return parcel_ptile_table_outliers(parcel_data, id_field, field, layer_name, inverse, nan_value)
    
    elif how == 'parcel_norm':
        return parcel_norm_table(parcel_data, id_field, field, layer_name, inverse, nan_value)
    
    else: #how = if_overlap
        return if_overlap (parcel_data=parcel_data, id_field=id_field, overlap_layer=join_data, field=field, layer_name=layer_name, points=points, overlap=overlap)
    



def category_merge(dfs: list, 
                   id_field, 
                   weights: dict, 
                   category_name:str):
    '''
    Pulls from a list of dfs and a "weights" dictionary to perform a weighted average of desired indicators.

    Inputs:
    dfs = list of data frames to include 
    weights = dictionary of variables and associated weights
    category name = a field name for the score variable 

    Outputs:
    A weighted sum (category name) and percentile ranking ('rnk_categoryname') for 
    each row based on given indicators and weights.

    Example:
    Define inputs
    local_accessibility_df = [schools, 
                            walkscore, 
                            towncenter]

    local_accessibility_weights = {
                                'nrm_schls': 1,     #school walksheds
                                'nrm_wlkscr': 2,    #walkscore
                                'nrm_twncntr': 2,   #within a town center
                                }

    Run function
    local_accessibility = category_merge(local_accessibility_df, local_accessibility_weights, 'lcl_acc')

    '''
    #merge all of the listed dataframes together based on their parcel id and geometry
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=[id_field, 'geometry'], how='outer'), dfs)

    #get weighted sum based on variable weights
    df_merged[category_name] = (df_merged[weights.keys()] * weights).sum(1) / sum(weights.values())

    #get percentile ranking for weighted sum value
    df_merged[('rnk_' + category_name)] = df_merged[category_name].rank(method='min', pct=True)

    #move geometry field to end
    df_merged.insert((len(df_merged.columns) - 1), 'geometry', df_merged.pop('geometry'))

    #move category score fields to beginning
    df_merged.insert(0, category_name, df_merged.pop(category_name))
    df_merged.insert(0, ('rnk_' + category_name), df_merged.pop('rnk_' + category_name))

    return(df_merged)



def comp_hists(col_1, col_2, muni, df, category):
    '''
    col_1 = raw data variable
    col_2 = procesed variable
    category = essentially a name for the title/exported chart
    '''
    #get IQR of col_1
    Q1=df[col_1].quantile(0.25)
    Q3=df[col_1].quantile(0.75)
    IQR=Q3-Q1

    val_1 = Q1-1.5*IQR
    val_2 = Q3+1.5*IQR
    
    sns.set_style("darkgrid")

    charts_fp = "K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Charts\\20230510_Histograms"
    charts_muni_fp = os.path.join(charts_fp, muni) 
    os.makedirs(charts_muni_fp, exist_ok=True) 

    plt.rcParams["figure.figsize"] = [10, 4]
    plt.rcParams["figure.autolayout"] = True

    fig, axes = plt.subplots(1, 3)

    fontsize = 9
    sns.histplot(df[col_1], ax=axes[0], bins=20, color='midnightblue').set_title('Raw value: Full distribution', fontsize=fontsize)
    sns.histplot(df[col_1], ax=axes[1], bins=20, color='midnightblue', binrange=(val_1, val_2)).set_title('Raw value: Normal data range (excluding outliers)', fontsize=fontsize)
    sns.histplot(df[col_2], ax=axes[2], bins=20, color='orange').set_title('Indicator value: Distribution after processing', fontsize=fontsize)
    fig.suptitle('Indicator: ' + category).set_fontweight('semibold')

    plt.savefig(charts_muni_fp + '\\' + category + '_hist.png')
    plt.close(fig)  

