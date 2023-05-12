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

def overlap_sjoin (target_layer, overlap_layer, field:str, stats=None):

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