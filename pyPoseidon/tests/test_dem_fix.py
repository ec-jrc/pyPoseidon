import pyPoseidon.dem as pdem
import pyPoseidon.grid as pg
import pytest
import numpy as np
import os
import geopandas as gp
import cartopy.feature as cf

from . import DATA_DIR

DEM_SOURCE = DATA_DIR / "dem.nc"


#define the lat/lon window and time frame of interest
window1 = {
    'lon_min' : -30,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 70.,
    'dem_source' : DEM_SOURCE,
}

window2={
     'lon_min':175., # lat/lon window
     'lon_max':185.,
     'lat_min':-21.5,
     'lat_max':-14.5,
     'dem_source' : DEM_SOURCE,
}

window3 = {
     'lon_min':-185., # lat/lon window
     'lon_max':-175.,
     'lat_min':-21.5,
     'lat_max':-14.5,
     'dem_source' : DEM_SOURCE,
}



cr = 'l'
coast = cf.NaturalEarthFeature(
    category='physical',
    name='land',
    scale='{}m'.format({'l':110, 'i':50, 'h':10}[cr]))


natural_earth = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])    

coast = cf.GSHHSFeature(
    scale='auto', #'coarse', 'low', 'intermediate', 'high, or 'full'
    levels = [1])

GSHHS = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])    



@pytest.mark.parametrize('dic', [ window1 , window2, window3])
def test_answer(tmpdir, dic):
    # Just elevation  
    df = pdem.dem(**dic) #get dem 
    df.adjust(natural_earth)
    
    c1 = 'adjusted' in df.Dataset.data_vars
    
    # Schism grid
    grid_file = DATA_DIR / 'hgrid.gr3'
    grid = pg.grid(type = 'tri2d',grid_file=grid_file) # read grid
    xg = grid.Dataset.SCHISM_hgrid_node_x.values
    yg = grid.Dataset.SCHISM_hgrid_node_y.values
    dic.update({'grid_x':xg, 'grid_y':yg})
    
    df = pdem.dem(**dic) #get dem 
    df.adjust(natural_earth)
       
    c2 = 'fval' in df.Dataset.data_vars
    
    ## D3D grid
    dic.update({'resolution':.1})
    grid = pg.grid(type = 'r2d', **dic) 
    gr = grid.Dataset
    xp,yp=gr.lons, gr.lats
     
    dic.update({'grid_x':xp, 'grid_y':yp})
    
    #get dem 
    df = pdem.dem(**dic)   
    df.adjust(natural_earth)
    
    c3 = 'fval' in df.Dataset.data_vars
    
    
    
    assert all([c1,c2,c3])
