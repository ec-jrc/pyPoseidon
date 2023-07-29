import pyposeidon.meteo as pmeteo
import pytest
import xarray as xr
import os
import glob
import numpy as np
import pandas as pd
import shutil


from . import DATA_DIR

filename = (DATA_DIR / "meteo.nc").as_posix()


@pytest.mark.parametrize("name", [filename])
def test_meteo(tmpdir, name):
    rpath = str(tmpdir) + "/"
    d = pmeteo.Meteo(filename)
    d.to_output(solver_name="schism", rpath=rpath, meteo_split_by="D")
    d.to_output(solver_name="schism", rpath=rpath, filename="all.nc")

    # read schism meteo files
    files = glob.glob(rpath + "sflux/*.nc")
    files.sort()
    ma = []
    for ifile in files:
        g = xr.open_dataset(ifile)
        ts = "-".join(g.time.attrs["base_date"].astype(str)[:3])
        time_r = pd.to_datetime(ts)
        times = time_r + pd.to_timedelta(g.time.values, unit="D").round("H")
        g = g.assign_coords({"time": times})
        ma.append(g)

    b = xr.merge(ma)
    b.close()

    tlist = pd.to_datetime(b.time.data) - pd.to_datetime(b.time.data[0])  # convert to Schism's time coords
    tlist = tlist / pd.to_timedelta("1D")

    b = b.assign_coords({"time": tlist})

    al = xr.open_dataset(rpath + "sflux/all.nc")

    assert b.equals(al)
