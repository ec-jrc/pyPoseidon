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
    df = pm.meteo(meteo_source=filenames, engine='cfgrib', combine_by='nested', combine_forecast=True, xr_kwargs = {'concat_dim' : 'step'}) # use combine
    df0 = pm.meteo(meteo_source=[filenames[0]],engine='cfgrib') # each one seperately
    df1 = pm.meteo(meteo_source=[filenames[1]],engine='cfgrib') # each one seperately
    df2 = pm.meteo(meteo_source=[filenames[2]],engine='cfgrib') # each one seperately
    df3 = pm.meteo(meteo_source=[filenames[3]],engine='cfgrib') # each one seperately

    #merge the single files
    joined = xr.concat([df0.Dataset.isel(time=slice(0,12)), df1.Dataset.isel(time=slice(0,12)), df2.Dataset.isel(time=slice(0,12)), df3.Dataset], dim='time')

    assert joined.equals(df.Dataset) # compare
