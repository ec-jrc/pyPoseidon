import pyPoseidon.meteo as pm
import pytest
import xarray as xr
import os
from glob import glob
import numpy as np


def cfgrib():

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
    filenames = glob(path + '/uvp_*')
    filenames.sort()
    
    #read meteo files
    df = pm.meteo(meteo_files=filenames, combine_by='nested', combine_forecast=True, xr_kwargs = {'concat_dim' : 'step'}) # use combine
    df0 = pm.meteo(meteo_files=[filenames[0]]) # each one seperately
    df1 = pm.meteo(meteo_files=[filenames[1]]) # each one seperately
    df2 = pm.meteo(meteo_files=[filenames[2]]) # each one seperately
    df3 = pm.meteo(meteo_files=[filenames[3]]) # each one seperately
    
    #merge the single files
    joined = xr.concat([df0.Dataset.isel(time=slice(0,12)), df1.Dataset.isel(time=slice(0,12)), df2.Dataset.isel(time=slice(0,12)), df3.Dataset], dim='time')
    
    return joined.equals(df.Dataset) # compare
        
    
def test_answer():
    assert cfgrib() == True
    
