import pytest
import pyposeidon
import os
import multiprocessing

from . import DATA_DIR


MESH_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()
METEO_FILE = (DATA_DIR / "erai.grib").as_posix()

OPTS = pytest.mark.parametrize(
    "flag, dem, meteo, mesh",
    [
        pytest.param(["all"], DEM_FILE, METEO_FILE, MESH_FILE, id="all modules"),
        pytest.param(["dem"], DEM_FILE, None, MESH_FILE, id="only bathymetry"),
        pytest.param(["meteo"], None, METEO_FILE, None, id="only meteo"),
        pytest.param(["model"], None, None, None, id="only model setup"),
        pytest.param(["model", "meteo"], None, METEO_FILE, None, id="model + meteo"),
        pytest.param(["model", "dem"], DEM_FILE, None, MESH_FILE, id="model + dem"),
        pytest.param(["dem", "meteo"], DEM_FILE, METEO_FILE, MESH_FILE, id="dem + meteo"),
    ],
)


@OPTS
def test_update_selector(tmpdir, flag, dem, meteo, mesh):
    # define in a dictionary the properties of the model..
    dic = {
        "solver_name": "schism",
        "mesh_file": mesh,
        "manning": 0.12,
        "windrot": 0.00001,
        "tag": "test",
        "start_date": "2017-10-1 0:0:0",
        "time_frame": "12H",
        "meteo_source": meteo,  # meteo file
        "dem_source": dem,
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

    # initialize a model
    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model

    dic.update({"update": flag})

    b = pyposeidon.model.set(**dic)

    b.create()

    assert isinstance(b, pyposeidon.schism.Schism) == True
