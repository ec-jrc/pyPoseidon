import pyPoseidon.meteo as pm
import pytest
import xarray as xr
import os
from glob import glob
import numpy as np

from . import DATA_DIR


def test_answer():
    filenames = sorted(DATA_DIR.glob("uvp_*"))

    #read meteo files
    df = pm.meteo(meteo_source=filenames, meteo_engine='cfgrib', meteo_combine_by='nested', meteo_merge='last', meteo_xr_kwargs = {'concat_dim' : 'step'}) # use combine
    df0 = pm.meteo(meteo_source=[filenames[0]],meteo_engine='cfgrib') # each one seperately
    df1 = pm.meteo(meteo_source=[filenames[1]],meteo_engine='cfgrib') # each one seperately
    df2 = pm.meteo(meteo_source=[filenames[2]],meteo_engine='cfgrib') # each one seperately
    df3 = pm.meteo(meteo_source=[filenames[3]],meteo_engine='cfgrib') # each one seperately

    #merge the single files
    joined = xr.concat([df0.Dataset.isel(time=slice(0,12)), df1.Dataset.isel(time=slice(0,12)), df2.Dataset.isel(time=slice(0,12)), df3.Dataset], dim='time')


    df_ = pm.meteo(meteo_source=filenames, meteo_engine='cfgrib', meteo_combine_by='nested', meteo_merge='first', meteo_xr_kwargs = {'concat_dim' : 'step'}) # use combine

    #merge the single files
    joined_ = xr.concat([df0.Dataset.isel(time=slice(0,13)), df1.Dataset.isel(time=slice(1,13)), df2.Dataset.isel(time=slice(1,13)), df3.Dataset.isel(time=slice(1,None))], dim='time')

    assert joined.equals(df.Dataset) # compare
    assert joined_.equals(df_.Dataset) # compare
