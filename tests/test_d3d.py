import pytest
import pyposeidon
from pyposeidon.utils import data
import os
import multiprocessing

from . import DATA_DIR

# define in a dictionary the properties of the model..
case1 = {
    "geometry": {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
    "start_date": "2018-10-1",
    "time_frame": "12H",
    "solver_name": "d3d",
    "resolution": 0.2,  # grid resoltuion
    "map_step": 20,  # step for output of map field in d3d
    "restart_step": 60,  # when to output restart file
    "meteo_source": [(DATA_DIR / "uvp_2018100100.grib").as_posix()],
    "dem_source": (DATA_DIR / "dem.nc").as_posix(),
    #     'update':['all'] # optional to select update quantities
}


@pytest.mark.delft
@pytest.mark.parametrize("dic", [case1])
def test_d3d(tmpdir, dic):
    # initialize a model
    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model
    b = pyposeidon.model.set(**dic)
    b.execute()
    out_b = data.get_output(**dic)

    a = pyposeidon.read(rpath + "d3d_model.json")  # read model
    a.execute()
    out_a = data.get_output(**dic)

    assert out_a == out_b
