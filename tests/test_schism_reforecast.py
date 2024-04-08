import pytest
import pyposeidon
from pyposeidon.utils import cast, data
import pyposeidon.meteo as pm
import json
import pandas as pd
import datetime
import os
import numpy as np
import xarray as xr

from . import DATA_DIR


MESH_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()
METEO_FILES_1 = [(DATA_DIR / "uvp_2018100100.grib").as_posix()]
METEO_FILES_2 = [
    (DATA_DIR / name).as_posix()
    for name in (
        "uvp_2018100100.grib",
        "uvp_2018100112.grib",
        "uvp_2018100200.grib",
        "uvp_2018100212.grib",
    )
]


# define in a dictionary the properties of the model..
case = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "schism",
    "start_date": "2018-10-1 0:0:0",
    "time_frame": "12H",
    "dem_source": DEM_FILE,
    "meteo_source": METEO_FILES_1,
    "meteo_merge": "first",  # combine meteo
    "meteo_combine_by": "nested",
    "meteo_xr_kwargs": {"concat_dim": "step"},
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 0.5,
        "nhot": 1,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
    "scribes": 2,
}


# define in a dictionary the properties of the model..
check = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "schism",
    "start_date": "2018-10-1 0:0:0",
    "time_frame": "36H",
    "dem_source": DEM_FILE,
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 1.5,
        "nhot": 0,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
    "scribes": 2,
}


@pytest.mark.schism
def test_schism_reforecast_workflow(tmpdir):
    # initialize a model
    rpath = str(tmpdir) + "/schism/"
    case.update({"rpath": rpath + "20181001.00/"})  # use tmpdir for running the model

    b = pyposeidon.model.set(**case)

    b.execute()

    # creating a time sequence of the runs
    start_date = pd.to_datetime("2018-10-1 0:0:0")
    end_date = pd.to_datetime("2018-10-2 0:0:0")
    date_list = pd.date_range(start_date, end_date, freq="12H")

    # creating a sequence of folder to store the runs. In this case we name them after the date attribute.
    # NOTE that the first folder is the fisrt run already perfomed!!
    rpaths = [rpath + datetime.datetime.strftime(x, "%Y%m%d.%H") + "/" for x in date_list]

    # creating a sequence of folder from which we read the meteo.
    meteo = []
    for date in date_list:
        prev_date = pd.to_datetime(date) - pd.to_timedelta("12H")
        prev_date = prev_date.strftime(format="%Y-%m-%d %H:%M:%S")
        dr = pd.date_range(prev_date, date, freq="12H")
        names = ["uvp_" + datetime.datetime.strftime(x, "%Y%m%d%H") + ".grib" for x in dr]
        dur = [(DATA_DIR / name).as_posix() for name in names]
        meteo.append(dur)

    # set cast
    for l in range(len(rpaths) - 1):
        h = cast.set(
            solver_name="schism",
            model=b,
            ppath=rpaths[l],
            cpath=rpaths[l + 1],
            meteo=meteo[l + 1],
            sdate=date_list[l + 1],
        )
        h.run(execute=True)  # execute

    # Run check case - Total duration
    check.update({"rpath": rpath + "check/"})  # use tmpdir for running the model

    # Combine meteo appropriately

    m1 = pm.Meteo(meteo_source=METEO_FILES_2[0])
    m2 = pm.Meteo(meteo_source=METEO_FILES_2[1])
    m3 = pm.Meteo(meteo_source=METEO_FILES_2[2])
    m4 = pm.Meteo(meteo_source=METEO_FILES_2[3])

    # extract correct chunk

    w1 = m1.Dataset.isel(time=slice(0, 13))
    w2 = m2.Dataset.isel(time=slice(1, 13))  # note that we keep the 12 hour from the previous file
    w3 = m3.Dataset.isel(time=slice(1, 13))
    w4 = m4.Dataset.isel(time=slice(1, 13))

    # combine
    meteo = xr.combine_by_coords([w1, w2, w3, w4], combine_attrs="override")
    # saving
    check.update({"meteo_source": meteo})

    c = pyposeidon.model.set(**check)

    c.execute()

    # COMPARE
    output = data.get_output(folders=rpaths, solver_name="schism")

    total = data.get_output(folders=[rpath + "check/"], solver_name="schism")

    r = output.Dataset.isel(time=slice(0, 36))

    rb = []
    for var in total.Dataset.data_vars:
        if not total.Dataset[var].equals(r[var]):
            rb.append(var)

    print(rb)

    #    flag = True TODO
    #    for var in rb:
    #        flag = False
    #        mdif = np.abs(total.results.Dataset[var].values - output.results.Dataset[var].values).max()
    #        if mdif < 1.e-14 :
    #            flag = True
    #    print(mdif)
    expected = ["wetdry_side", "wetdry_elem", "wetdry_node", "zcor", "elev", "hvel"]
    assert rb == expected
