import pytest
import pyposeidon
import os
import multiprocessing
import geopandas as gp

from . import DATA_DIR

DEM_FILE = DATA_DIR / "dem.nc"

COAST_FILE = (DATA_DIR / "ocean.zip").as_posix()

WINDOWS = pytest.mark.parametrize(
    "window",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 165.0, "lon_max": 185.0, "lat_min": -50, "lat_max": -30},
        {"lon_min": -185.0, "lon_max": -165.0, "lat_min": -50, "lat_max": -30},
        {"lon_min": -25.0, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 68.0},
    ],
)


@pytest.mark.schism
@WINDOWS
def test_schism(tmpdir, window):
    # initialize a model
    dic = {
        "solver_name": "schism",
        "geometry": window,
        "coastlines": COAST_FILE,
        "manning": 0.12,
        "windrot": 0.00001,
        "tag": "test",
        "start_date": "2017-10-1 0:0:0",
        "time_frame": "12h",
        "mesh_generator": "gmsh",
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
        "scribes": 1,
    }

    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model

    b = pyposeidon.model.set(**dic)

    b.execute()
    b.results()

    err_file = b.rpath + "/outputs/fatal.error"
    assert os.stat(err_file).st_size == 0
