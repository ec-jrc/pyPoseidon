import pytest
import pyposeidon
import os
import subprocess
import pyposeidon.model as pm
import pyposeidon.mesh as pmesh

from . import DATA_DIR


MESH_FILE = (DATA_DIR / "hgrid.gr3").as_posix()
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

# define in a dictionary the properties of the model..
case1 = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "start_date": "2017-10-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "erai.grib").as_posix()],  # meteo file
    "update": ["model"],  # update only model
}

case2 = {
    "solver_name": "schism",
    "mesh_file": MESH_FILE,
    "epath": "/wrong_path/",
    "start_date": "2017-10-1 0:0:0",
    "time_frame": "12H",
    "meteo_source": [(DATA_DIR / "erai.grib").as_posix()],  # meteo file
    "update": ["model"],  # update only model
}

geometry = {
    "lon_min": -25.0,  # lat/lon window
    "lon_max": -9.0,
    "lat_min": 56.0,
    "lat_max": 74.0,
}


def jigsaw(tmpdir):
    rpath = str(tmpdir) + "/"
    mesh = pmesh.set(
        type="tri2d",
        geometry=geometry,
        rpath=rpath,
        mesh_generator="jigsaw",
        hfun_max=-1,
    )


def gmsh(tmpdir):
    rpath = str(tmpdir) + "/"
    mesh = pmesh.set(
        type="tri2d",
        geometry=geometry,
        rpath=rpath,
        mesh_generator="gmsh",
        use_bindings=False,
        gmsh_args={"clmax": -2},
    )


@pytest.mark.parametrize("case", [case1, case2])
def test_schism_execution(tmpdir, case):
    # initialize a model
    rpath = str(tmpdir) + "/"
    case.update({"rpath": rpath})  # use tmpdir for running the model
    model = pm.set(**case)
    with pytest.raises(subprocess.CalledProcessError) as exc:
        model.execute()


def test_jigsaw(tmpdir):
    with pytest.raises(subprocess.CalledProcessError) as exc:
        jigsaw(tmpdir)
    assert "error" in str(exc)


def test_gmsh(tmpdir):
    with pytest.raises(subprocess.CalledProcessError) as exc:
        gmsh(tmpdir)
    assert "Error" in str(exc)
