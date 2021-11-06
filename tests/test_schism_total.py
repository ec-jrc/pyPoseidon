import pytest
import pyposeidon
import os
import multiprocessing

from . import DATA_DIR

DEM_FILE = DATA_DIR / "dem.nc"

case0 = {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0}

case1 = {
    "lon_min": 175.0,
    "lon_max": 184.0,
    "lat_min": -21.5,
    "lat_max": -14.5,
}  # lat/lon window

case2 = {"lon_min": -20.0, "lon_max": -10.0, "lat_min": 63.0, "lat_max": 67.0}

case3 = {
    "lon_min": 175.0 - 360.0,
    "lon_max": 184.0 - 360.0,
    "lat_min": -21.5,
    "lat_max": -14.5,
}  # lat/lon window

case4 = {"lon_min": -25.0, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 68.0}


@pytest.mark.schism
@pytest.mark.parametrize("case", [case0, case2])  # , case1, case3])
def test_schism(tmpdir, case):
    # initialize a model
    dic = {
        "solver": "schism",
        "geometry": case,
        "manning": 0.12,
        "windrot": 0.00001,
        "tag": "test",
        "start_date": "2017-10-1 0:0:0",
        "time_frame": "12H",
        "grid_generator": "jigsaw",
        "meteo_source": [DATA_DIR / "erai.grib"],  # meteo file
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

    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model

    b = pyposeidon.model.set(**dic)

    b.execute()
    b.results()

    err_file = b.rpath + "/outputs/fatal.error"
    assert os.stat(err_file).st_size == 0
