import pyPoseidon.grid as pg
import numpy as np
import pytest
import os
import geopandas as gp
import cartopy.feature as cf


cr='l'
coast = cf.NaturalEarthFeature(
    category='physical',
    name='land',
    scale='{}m'.format({'l':110, 'i':50, 'h':10}[cr]))


natural_earth = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])    

coast = cf.GSHHSFeature(
    scale='auto',
    levels = [1])

GSHHS = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])    




#define the lat/lon window and time frame of interest
window0 = {
    'lon_min' : -30,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 70.}

window1={'lon_min':175., # lat/lon window
     'lon_max':184.,
     'lat_min':-21.5,
     'lat_max':-14.5}
     
window2 = {
     'lon_min' : -20.,
     'lon_max' : -10.,
     'lat_min' : 63.,
     'lat_max' : 67.}

window3={'lon_min':175.-360., # lat/lon window
     'lon_max':184.-360.,
     'lat_min':-21.5,
     'lat_max':-14.5}
     
window4 = {
    'lon_min' : -25.,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 68.}

@pytest.mark.parametrize('window', [ window0, window1, window2, window3, window4 ])
@pytest.mark.parametrize('coast', [ natural_earth, GSHHS ])

def test_answer(tmpdir, window, coast):
    
    df = pg.grid(type='tri2d',geometry=window, coastlines=coast, rpath = str(tmpdir)+'/')
    
    check = np.isnan(df.Dataset.depth.values).sum() == 0
    assert check == True