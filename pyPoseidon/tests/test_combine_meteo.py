import pyPoseidon.meteo as pmeteo
import pytest
import xarray as xr
import os
from glob import glob
import numpy as np


def cfgrib():

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
    filenames = glob(path + '/oper/*')
    filenames.sort()
    
    #read meteo files
    df = pmeteo.cfgrib(filenames=filenames[1:], combine=True) # use combine
    df1 = pmeteo.cfgrib(filenames=filenames[1]) # each one seperately
    df2 = pmeteo.cfgrib(filenames=filenames[2]) # each one seperately
    
    #merge the single files
    joined = xr.concat([df1.isel(time=slice(0,12)), df2], dim='time')
    
    return joined.equals(df) # compare
        
    
def test_answer():
    assert cfgrib() == True
    
