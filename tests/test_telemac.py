import os
import pathlib

import pyposeidon
import pytest

from . import DATA_DIR

try:
    import telapy
except ImportError:
    pytest.skip(reason="requires telapy", allow_module_level=True)


MESH_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

# define in a dictionary the properties of the model..
case0 = {
    "solver_name": "telemac",
    "mesh_file": MESH_FILE,
    "start_date": "2017-10-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "erai.grib").as_posix()],  # meteo file
    "dem_source": DEM_FILE,
    "update": ["all"],  # update only meteo, keep dem
    "dt": 50,
    "parameters": {
        "chezy": 30,
    },
}

case1 = {
    "solver_name": "telemac",
    "mesh_file": MESH_FILE,
    "start_date": "2017-10-1 0:0:0",
    "end_date": "2017-10-2 0:0:0",  # optional instead of time_frame
    "dem_source": DEM_FILE,
    "meteo_source": [(DATA_DIR / "erai.grib").as_posix()],  # meteo file
    "update": ["all"],  # update only meteo, keep dem
    "dt": 100,
    "parameters": {
        "chezy": 100,
    },
}


case2 = {
    "solver_name": "telemac",
    "mesh_file": MESH_FILE,
    "start_date": "2011-1-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "era5.grib").as_posix()],  # meteo file
    "dem_source": DEM_FILE,
    "update": ["all"],  # update only meteo, keep dem
    "dt": 200,
}


case3 = {
    "solver_name": "telemac",
    "mesh_file": MESH_FILE,
    "module": "tomawac",
    "start_date": "2011-1-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "era5.grib").as_posix()],  # meteo file
    "dem_source": DEM_FILE,
    "update": ["all"],  # update only meteo, keep dem
    "dt": 1800,
}


def telemac(tmpdir, dic):
    # initialize a model
    rpath = str(tmpdir) + "/"
    dic.update({"rpath": rpath})  # use tmpdir for running the model

    b = pyposeidon.model.set(**dic)

    b.execute()
    b.results()
    zarr_tar = b.rpath + "/outputs/out_2D.zarr.tar"
    if os.path.exists(zarr_tar):
        a = pyposeidon.model.read(rpath + b.module + "_model.json")  # read model
        a.execute()
        a.results()
        zarr_tar = a.rpath + "/outputs/out_2D.zarr.tar"
        return os.path.exists(zarr_tar)
    else:
        return False


@pytest.mark.telemac
@pytest.mark.parametrize("case", [case0, case1, case2, case3])
def test_answer(tmpdir, case):
    assert telemac(tmpdir, case) == True
