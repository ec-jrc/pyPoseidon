import pytest
import pyposeidon
import os
import multiprocessing
import geopandas as gp
import cartopy.feature as cf

from . import DATA_DIR

DEM_FILE = DATA_DIR / "dem.nc"

coast = cf.GSHHSFeature(scale="auto", levels=[1])

GSHHS = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])


WINDOWS = pytest.mark.parametrize(
    "window",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 175.0, "lon_max": 185.0, "lat_min": -21.5, "lat_max": -14.5},
        {"lon_min": -185.0, "lon_max": -175.0, "lat_min": -21.5, "lat_max": -14.5},
        {"lon_min": -25.0, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 68.0},
    ],
)


@pytest.mark.schism
@WINDOWS
def test_schism(tmpdir, window):
    # initialize a model
    dic = {
        "solver": "schism",
        "geometry": window,
        "coastlines": GSHHS,
        "manning": 0.12,
        "windrot": 0.00001,
        "tag": "test",
        "start_date": "2017-10-1 0:0:0",
        "time_frame": "12H",
        "mesh_generator": "jigsaw",
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
