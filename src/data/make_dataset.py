# -*- coding: utf-8 -*-
import click
import logging
import geopandas as gpd
from pathlib import Path
import pandas as pd
import os
#from dotenv import find_dotenv, load_dotenv

def get_landuse_data(muni):
    '''
    input = muni name
    output = picks out the right shapefile from the state's municipal land use database;
             makes a subdirectory in intermediate folder w town name and exports land use shapefile to it
             reads that shapefile in as a geodataframe

    '''
    intermediate_path = "K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Data\\Intermediate"
    mapc_lpd_fp = "K:\\DataServices\\Projects\\Current_Projects\\Housing\\Section_3A\\Analytical_Toolbox\\Project_Files\\Parcel_Database\\MAPC_LandParcelDatabase.csv"
    mapc_lpd = pd.read_csv(mapc_lpd_fp)

    path = os.path.join(intermediate_path, muni)


    #set land  use variable
    town_lu_path = ""
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if filename.endswith('.shp'):
                town_lu_path = os.path.join(dirpath, filename)

    muni_state_parcels = gpd.read_file(town_lu_path).rename(columns={'LOC_ID': 'parloc_id'})
    
    muni_lpd_preprocess = muni_state_parcels.merge(mapc_lpd, on='parloc_id', how='left')

    #unconstrained land area

    #create a non-excluded land field, convert to acres 
    muni_lpd_preprocess['non_exc'] = (muni_lpd_preprocess['SQFT'] - muni_lpd_preprocess['Tot_Exclud']) / 43560
    
    #total amount of excluded acres
    muni_lpd_preprocess['exc_acres'] = muni_lpd_preprocess['Tot_Exclud'] / 43560

    #total land are in acres
    muni_lpd_preprocess['acreage'] = muni_lpd_preprocess['SQFT'] / 43560

    #total amount of developable land
    muni_lpd_preprocess['developable_land'] = muni_lpd_preprocess['acreage'] - muni_lpd_preprocess['exc_acres']
    
    #percent of land area that is not excluded
    muni_lpd_preprocess['pct_developable'] = 1 - (muni_lpd_preprocess['exc_acres'] / muni_lpd_preprocess['acreage'])

    #muni_lpd_preprocess = muni_lpd_preprocess.rename(columns={'LOC_ID': 'parloc_id'})

    return(muni_lpd_preprocess)  
    