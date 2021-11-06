import pytest
import pyposeidon
import os
import multiprocessing

from . import DATA_DIR


GRID_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

# define in a dictionary the properties of the model..
case1 = {
    "solver": "schism",
    "grid_file": GRID_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "test",
    "start_date": "2017-10-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "erai.grib").as_posix()],  # meteo file
    "dem_source": DEM_FILE,
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 0.3,
        "nhot": 0,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
}

case2 = {
    "solver": "schism",
    "grid_file": GRID_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "test",
    "start_date": "2017-10-1 0:0:0",
    "end_date": "2017-10-2 0:0:0",  # optional instead of time_frame
    "dem_source": DEM_FILE,
    "meteo_source": [(DATA_DIR / "erai.grib").as_posix()],  # meteo file
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 0.3,
        "nhot": 0,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
}


case3 = {
    "solver": "schism",
    "grid_file": GRID_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "test",
    "start_date": "2011-1-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "era5.grib").as_posix()],  # meteo file
    "dem_source": DEM_FILE,
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 0.3,
        "nhot": 0,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
}


def schism(tmpdir, dic):
    # initialize a model
    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model

    b = pyposeidon.model.set(**dic)

    b.execute()
    b.results()
    err_file = b.rpath + "/outputs/fatal.error"
    if os.stat(err_file).st_size == 0:
        a = pyposeidon.model.read(rpath + "test_model.json")  # read model
        a.execute()
        a.results()
        err_file = a.rpath + "/outputs/fatal.error"
        return os.stat(err_file).st_size == 0
    else:
        return False


@pytest.mark.schism
@pytest.mark.parametrize("case", [case1, case2, case3])
def test_answer(tmpdir, case):
    assert schism(tmpdir, case) == True
