import pyposeidon.meteo as pmeteo
import pytest
import xarray as xr
import os
import numpy as np
import shutil


from . import DATA_DIR

filename = (DATA_DIR / 'meteo.nc').as_posix()

@pytest.mark.parametrize('name', [filename])
def test_meteo(name):
    try:
        d = pmeteo.meteo(filename, meteo_engine='netcdf')
        r = True
    except:
        r = False

    assert r == True
