import os
import numpy as np
import xarray as xr
import multiprocessing

# from . import DATA_DIR
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


@pytest.mark.telemac
@pytest.mark.parametrize("case", [case0, case1, case2, case3])
def test_telemac(tmpdir, case):
    # initialize a model
    rpath = str(tmpdir) + "/"
    case.update({"rpath": rpath})  # use tmpdir for running the model

    b = pyposeidon.model.set(**case)
    b.create()
    b.mesh.Dataset.type[:] = "closed"
    b.output()
    b.save()
    b.run()
    b.results()
    zarr_tar = b.rpath + "/outputs/out_2D.zarr.tar"
    assert os.path.exists(zarr_tar)
    if os.path.exists(zarr_tar):
        a = pyposeidon.model.read(rpath + b.module + "_model.json")  # read model
        a.create()
        a.mesh.Dataset.type[:] = "closed"
        a.output()
        a.save()
        a.run()
        a.results()
        zarr_tar = a.rpath + "/outputs/out_2D.zarr.tar"
        res = a.rpath + "/results_2D.slf"
        assert os.path.exists(zarr_tar)
        ds = xr.open_dataset(res, engine="selafin")
        if "WH" in ds.variables:
            max_WH = ds.WH.max().values
            assert np.isclose(max_WH, 4.13, atol=1e-2)
        elif "S" in ds.variables:
            max_S = ds.S.max().values
            print(max_S)
            assert max_S < 1
        else:
            pass
