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
    df = pmeteo(meteo_files=filenames, combine=True) # use combine
    df0 = pmeteo(meteo_files=[filenames[0]]) # each one seperately
    df1 = pmeteo(meteo_files=[filenames[1]]) # each one seperately
    df2 = pmeteo(meteo_files=[filenames[2]]) # each one seperately
    df3 = pmeteo(meteo_files=[filenames[3]]) # each one seperately
    
    #merge the single files
    joined = xr.concat([df0.uvp.isel(time=slice(0,12)), df1.uvp.isel(time=slice(0,12)), df2.uvp.isel(time=slice(0,12)), df3.uvp], dim='time')
    
    return joined.equals(df.uvp) # compare
        
    
def test_answer():
    assert cfgrib() == True
    
