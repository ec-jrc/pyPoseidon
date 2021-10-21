import pytest
import pyposeidon
from pyposeidon.utils import data
import os
import multiprocessing

NCORES = max(1, multiprocessing.cpu_count() - 1)

from . import DATA_DIR

# define in a dictionary the properties of the model..
case1 = {
    "geometry": {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
    "start_date": "2018-10-1",
    "time_frame": "12H",
    "solver": "d3d",
    "resolution": 0.2,  # grid resoltuion
    "map_step": 20,  # step for output of map field in d3d
    "restart_step": 60,  # when to output restart file
    "ncores": NCORES,  # number of cores
    "meteo_source": [(DATA_DIR / "uvp_2018100100.grib").as_posix()],
    "dem_source": (DATA_DIR / "dem.nc").as_posix(),
    #     'update':['all'] # optional to select update quantities
}


def d3d(tmpdir, dic):
    # initialize a model
    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model
    b = pyposeidon.model(**dic)

    try:
        b.execute()
        out = data.data(**dic)
        a = pyposeidon.read_model(rpath + "d3d_model.json")  # read model
        a.execute()
        out = data.data(**dic)
        return True
    except:
        return False


@pytest.mark.delft
@pytest.mark.parametrize("case", [case1])
def test_answer(tmpdir, case):
    assert d3d(tmpdir, case) == True
