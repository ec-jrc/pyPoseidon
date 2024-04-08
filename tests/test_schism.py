import copy
import glob
import os

import pandas as pd
import pytest
import pyposeidon
import pyposeidon.schism

from . import DATA_DIR


MESH_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

# define in a dictionary the properties of the model..
case1 = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
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
    "scribes": 2,
}

case2 = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
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
    "scribes": 2,
}


case3 = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
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
    "scribes": 2,
}


case4 = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "test",
    "start_date": "2011-1-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "era5.grib").as_posix()],  # meteo file
    "dem_source": DEM_FILE,
    "monitor": True,
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 0.3,
        "nhot": 0,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
        "nc_out": 0,
    },
    "scribes": 2,
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
@pytest.mark.parametrize("case", [case1, case2, case3, case4])
def test_answer(tmpdir, case):
    assert schism(tmpdir, case) == True


def test_schism_meteo_split_by(tmpdir):
    model_description = copy.deepcopy(case1)
    model_description.update(
        {
            "rpath": tmpdir,
            "meteo_split_by": "1D",
            "time_frame": "144H",
            "update": ["meteo"],
        }
    )
    model = pyposeidon.model.set(**model_description)
    model.create()
    model.output()
    assert len(glob.glob(f"{tmpdir}/sflux/*.nc")) == 7, glob.glob(f"{tmpdir}/sflux/*.nc")


def test_parse_mirror_out():
    path = DATA_DIR / "mirror.out"
    df = pyposeidon.schism.parse_mirror_out(path)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 65
    assert "etatot" in df
    assert "etaavg" in df
    assert max(df.etatot) == pytest.approx(26.501, abs=1e-3)
    assert max(df.etaavg) == pytest.approx(0.0227, abs=1e-3)
    assert df.index[0] == pd.Timestamp(2017, 10, 1)
    assert df.index[-1] == pd.Timestamp(2017, 10, 1, 7, 6, 40)
