import pyPoseidon.dem as pdem
import pyPoseidon.grid as pg
import pytest
import numpy as np
import os
import geopandas as gp
import cartopy.feature as cf

from . import DATA_DIR

DEM_SOURCE = DATA_DIR / 'dem.nc'


#define the lat/lon window and time frame of interest
window1 = {
    'lon_min' : -30,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 70.
}

window2={
     'lon_min':175., # lat/lon window
     'lon_max':185.,
     'lat_min':-21.5,
     'lat_max':-14.5
}

window3 = {
     'lon_min':-185., # lat/lon window
     'lon_max':-175.,
     'lat_min':-21.5,
     'lat_max':-14.5
}

window4 = {
    'lon_min' : -3,
    'lon_max' : 1.,
    'lat_min' : 42.,
    'lat_max' : 45.
}



cr = 'i'
coast = cf.NaturalEarthFeature(
    category='physical',
    name='land',
    scale='{}m'.format({'l':110, 'i':50, 'h':10}[cr]))


natural_earth = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])

coast = cf.GSHHSFeature(
    scale='auto', #'coarse', 'low', 'intermediate', 'high, or 'full'
    levels = [1])

GSHHS = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])


@pytest.mark.parametrize('dic', [ window1 , window2, window3, window4])
def test_elevation(tmpdir, dic):
    # Just elevation
    df = pdem.dem(**dic, dem_source = DEM_SOURCE) #get dem
    df.adjust(natural_earth)

    assert np.isnan(df.Dataset.adjusted.values).sum() == 0

    
# Schism grid    
#@pytest.mark.schism
@pytest.mark.parametrize('dic', [ window1 , window2, window3, window4])
def test_schism_grid(tmpdir, dic):
    
    grid = pg.grid(type = 'tri2d',geometry=dic, coastlines=natural_earth, rpath = str(tmpdir)+'/')
    xg = grid.Dataset.SCHISM_hgrid_node_x.values
    yg = grid.Dataset.SCHISM_hgrid_node_y.values

    df = pdem.dem(**dic, grid_x=xg, grid_y=yg, dem_source=DEM_SOURCE) #get dem
    df.adjust(natural_earth) 
          
    assert np.isnan(df.Dataset.fval.values).sum() == 0
    
# D3D grid
@pytest.mark.delft
@pytest.mark.parametrize('dic', [ window1 , window2, window3, window4])
def test_d3d_grid(tmpdir, dic):
    
    grid = pg.grid(type = 'r2d', geometry=dic, resolution=.1, rpath = str(tmpdir)+'/')
    gr = grid.Dataset
    xp,yp=gr.lons, gr.lats

    #get dem
    df = pdem.dem(**dic, grid_x=xp, grid_y=yp, dem_source=DEM_SOURCE)
    df.adjust(natural_earth)

    assert np.isnan(df.Dataset.fval.values).sum() == 0
