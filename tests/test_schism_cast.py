from copy import deepcopy
import datetime
import filecmp
import os

import pytest
import pyposeidon
from pyposeidon.utils import cast, data
import pandas as pd

from . import DATA_DIR


MESH_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()
METEO_FILES_1 = [(DATA_DIR / name).as_posix() for name in ("uvp_2018100100.grib", "uvp_2018100112.grib")]
METEO_FILES_2 = [
    (DATA_DIR / name).as_posix()
    for name in (
        "uvp_2018100100.grib",
        "uvp_2018100112.grib",
        "uvp_2018100200.grib",
        "uvp_2018100212.grib",
    )
]


# define in a dictionary the properties of the model..
test_case = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "schism",
    "start_date": "2018-10-1 0:0:0",
    "time_frame": "12h",
    "dem_source": DEM_FILE,
    "meteo_source": METEO_FILES_1,
    "meteo_merge": "last",  # combine meteo
    "meteo_combine_by": "nested",
    "meteo_xr_kwargs": {"concat_dim": "step"},
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 0.5,
        "nhot": 1,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
    "scribes": 1,
}

# define in a dictionary the properties of the model..
check = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "manning": 0.12,
    "windrot": 0.00001,
    "tag": "schism",
    "start_date": "2018-10-1 0:0:0",
    "time_frame": "36h",
    "dem_source": DEM_FILE,
    "meteo_source": METEO_FILES_2,
    "meteo_merge": "last",  # combine meteo
    "meteo_combine_by": "nested",
    "meteo_xr_kwargs": {"concat_dim": "step"},
    "update": ["all"],  # update only meteo, keep dem
    "parameters": {
        "dt": 400,
        "rnday": 1.5,
        "nhot": 0,
        "ihot": 0,
        "nspool": 9,
        "ihfskip": 36,
        "nhot_write": 108,
    },
    "scribes": 1,
}


@pytest.mark.schism
@pytest.mark.parametrize(
    "copy",
    [pytest.param(True, id="copy files"), pytest.param(False, id="symlink files")],
)
def test_schism_cast(tmpdir, copy):
    base_rpath = os.path.join(tmpdir, "schism")
    original_rpath = os.path.join(base_rpath, "20181001.00")
    next_rpath = os.path.join(base_rpath, "20181001.12")

    # initialize a model
    original_model_info = deepcopy(test_case)
    original_model_info.update(
        {
            "rpath": original_rpath,
            "update": ["meteo"],
        },
    )
    original_model = pyposeidon.model.set(**original_model_info)
    original_model.create()  # constructs all required parts e.g. mesh, dem, meteo, etc.
    original_model.output()  # save to files
    original_model.save()  # saves the json model reference file
    original_model.set_obs()  # setup station points
    original_model.execute()

    casted = cast.set(
        solver_name="schism",
        model=original_model,
        ppath=original_model.rpath,
        cpath=next_rpath,
        meteo=(DATA_DIR / "uvp_2018100112.grib").as_posix(),
        sdate=pd.to_datetime("2018-10-01 12:00:00"),
        copy=copy,
    )
    casted.run(execute=False)

    copied_files = cast.SchismCast.files + cast.SchismCast.station_files
    if copy:
        copied_files.extend(cast.SchismCast.model_files)

    for filename in copied_files:
        original_file = os.path.join(original_model.rpath, filename)
        if os.path.exists(original_file):
            copied_file = os.path.join(next_rpath, filename)
            assert os.path.exists(copied_file)
            assert not os.path.islink(copied_file)
            assert filecmp.cmp(original_file, copied_file, shallow=False)

    if not copy:
        for filename in cast.SchismCast.model_files:
            original_file = os.path.join(original_model.rpath, filename)
            if os.path.exists(original_file):
                symlinked_file = os.path.join(next_rpath, filename)
                assert os.path.exists(symlinked_file)
                assert os.path.islink(symlinked_file)
                if not os.path.islink(original_file):
                    assert os.path.realpath(symlinked_file) == original_file


@pytest.mark.schism
@pytest.mark.parametrize(
    "hotstart",
    [1, 2],
)
def test_schism_cast_workflow(tmpdir, hotstart):
    # initialize a model
    rpath = str(tmpdir) + "/schism/"
    test_case.update({"rpath": rpath + "20181001.00/"})  # use tmpdir for running the model

    b = pyposeidon.model.set(**test_case)

    b.execute()

    # creating a time sequence of the runs
    start_date = pd.to_datetime("2018-10-1 0:0:0")
    end_date = pd.to_datetime("2018-10-2 0:0:0")
    date_list = pd.date_range(start_date, end_date, freq="12h")

    # creating a sequence of folder to store the runs. In this case we name them after the date attribute.
    # NOTE that the first folder is the fisrt run already perfomed!!
    rpaths = [rpath + datetime.datetime.strftime(x, "%Y%m%d.%H") + "/" for x in date_list]

    # creating a sequence of folder from which we read the meteo.
    meteo = []
    for date in date_list:
        end_date = pd.to_datetime(date) + pd.to_timedelta("12h")
        end_date = end_date.strftime(format="%Y-%m-%d %H:%M:%S")
        dr = pd.date_range(date, end_date, freq="12h")
        names = ["uvp_" + datetime.datetime.strftime(x, "%Y%m%d%H") + ".grib" for x in dr]
        dur = [(DATA_DIR / name).as_posix() for name in names]
        meteo.append(dur)

    # set cast
    for l in range(len(rpaths) - 1):
        h = cast.set(
            solver_name="schism",
            model=b,
            ppath=rpaths[l],
            cpath=rpaths[l + 1],
            meteo=meteo[l + 1],
            sdate=date_list[l + 1],
            ihot=hotstart,
        )

        h.run(execute=True)  # execute

    # Run check case - Total duration
    check.update({"rpath": rpath + "check/"})  # use tmpdir for running the model

    c = pyposeidon.model.set(**check)

    c.execute()

    # COMPARE
    output = data.get_output(folders=rpaths, solver_name="schism")

    total = data.get_output(folders=[rpath + "check/"], solver_name="schism")

    rb = []
    for var in total.Dataset.data_vars:
        if not total.Dataset[var].equals(output.Dataset[var]):
            rb.append(var)

    #    print(rb)

    #    flag = True TODO
    #    for var in rb:
    #        flag = False
    #        mdif = np.abs(total.results.Dataset[var].values - output.results.Dataset[var].values).max()
    #        if mdif < 1.e-14 :
    #            flag = True
    #    print(mdif)

    #    expected = ["wetdry_side", "wetdry_elem", "wetdry_node", "zcor", "elev", "hvel"]
    assert (rb == ["zcor"]) or (rb == [])
