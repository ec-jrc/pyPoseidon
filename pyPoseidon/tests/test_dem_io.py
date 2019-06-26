import pyPoseidon.dem as pdem
import pyPoseidon.model as pmodel
import pyPoseidon.grid as pgrid


import pytest
import xarray as xr
import os
import shutil
import numpy as np


#define the lat/lon window and time frame of interest
window1 = {
    'minlon' : -30,
    'maxlon' : -10.,
    'minlat' : 60.,
    'maxlat' : 70.
}



def schism(tmpdir,kwargs):

    grid = pgrid.grid(type='tri2d',grid_file='./data/hgrid.gr3')
    
    #update dem
    
    xp = grid.impl.Dataset.SCHISM_hgrid_node_x.values
    yp = grid.impl.Dataset.SCHISM_hgrid_node_y.values     
    
    kwargs.update({'grid_x':xp, 'grid_y':yp})
    #get dem 
    df = pdem.dem(**kwargs)

    new = df.altimetry.ival # get the interpolated values
    new = new.rename({'k':'nSCHISM_hgrid_node'}) # rename coord
    new.name = 'depth' # set correct name
    
    grid.impl.Dataset['depth'] = new # assign

    filename_ = str(tmpdir.join('hgrid_.gr3'))
    #output to grid file
    grid.impl.to_file(filename_)

    #read again new grid
    grid_ = pgrid.grid(type='tri2d',grid_file=filename_)

    #compare
    return grid.impl.Dataset.equals(grid_.impl.Dataset)
    
def d3d(tmpdir, kwargs):
    
    ## lat,lon grid
    resolution=.1
    lon=np.arange(kwargs['minlon'],kwargs['maxlon'],resolution)
    lat=np.arange(kwargs['minlat'],kwargs['maxlat'],resolution)
    xp, yp = np.meshgrid(lon,lat)
    
    kwargs.update({'grid_x':xp, 'grid_y':yp})
    
    #get dem 
    df = pdem.dem(**kwargs)

    rpath = str(tmpdir)+'/'
    #output 
    pdem.to_output(df.altimetry,solver='d3d',rpath=rpath)

    #read again dem
    rd = pmodel.d3d.from_dep(rpath + 'd3d.dep')

    #compare  
    c1 = -rd.where(rd!=-999)
    c2 = df.altimetry.ival.where(df.altimetry.ival < 0)

    return c1.fillna(0).equals(c2.fillna(0))
    
    
    
@pytest.mark.parametrize('kwargs', [ window1 ])
def test_answer(tmpdir, kwargs):
    
    
    assert schism(tmpdir,kwargs) == True
    assert d3d(tmpdir,kwargs) == True
    
