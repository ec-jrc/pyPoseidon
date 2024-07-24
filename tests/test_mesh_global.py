import os

import geopandas as gp
import numpy as np
import pyposeidon.mesh as pmesh
import pytest

from . import DATA_DIR

DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

COAST_FILE = (DATA_DIR / "ocean.zip").as_posix()


@pytest.mark.parametrize("ggor", ["jigsaw", "gmsh"])
@pytest.mark.parametrize("bgmesh", [None, DEM_FILE])
@pytest.mark.parametrize("bindings", [True, False])
@pytest.mark.parametrize("cbuffer", [None, 0.001])
def test_io(pytestconfig, tmpdir, ggor, bgmesh, bindings, cbuffer):
    # Skip the test unless --runslow has been passed
    if bgmesh is not None:
        if not pytestconfig.getoption("--runslow"):
            pytest.skip("slow test")

    mesh = pmesh.set(
        type="tri2d",
        geometry="global",
        coastlines=COAST_FILE,
        rpath=str(tmpdir) + "/",
        mesh_generator=ggor,
        dem_source=bgmesh,
        use_bindings=bindings,
        cbuffer=cbuffer,
    )

    # save to file
    filename = str(tmpdir.join("hgrid_.gr3"))
    mesh.to_file(filename)

    # read from file
    m = pmesh.set(type="tri2d", mesh_file=filename)

    dic = {}
    for d in m.Dataset.data_vars:
        dic.update({d: m.Dataset[d].equals(mesh.Dataset[d])})

    dic.pop("SCHISM_hgrid_node_x", None)
    dic.pop("SCHISM_hgrid_node_y", None)

    check1 = all(value == True for value in dic.values())

    dx = np.abs(m.Dataset.SCHISM_hgrid_node_x - mesh.Dataset.SCHISM_hgrid_node_x).max()
    dy = np.abs(m.Dataset.SCHISM_hgrid_node_y - mesh.Dataset.SCHISM_hgrid_node_y).max()

    check2 = (dx.values.max() < 1.0e-8) & (dy.values.max() < 1.0e-8)

    assert all([c == True for c in [check1, check2]])


@pytest.mark.schism
@pytest.mark.parametrize("ggor", ["jigsaw", "gmsh"])
@pytest.mark.parametrize("bgmesh", [None, DEM_FILE])
@pytest.mark.parametrize("bindings", [True, False])
@pytest.mark.parametrize("cbuffer", [None, 0.001])
def test_validate(pytestconfig, tmpdir, ggor, cbuffer, bgmesh, bindings):
    if bgmesh is not None:
        if not pytestconfig.getoption("--runslow"):
            pytest.skip("slow test")
    if ggor == "jigsaw":
        pytest.xfail("Fixing these is a WIP")

    mesh = pmesh.set(
        type="tri2d",
        geometry="global",
        coastlines=COAST_FILE,
        rpath=str(tmpdir) + "/",
        mesh_generator=ggor,
        dem_source=bgmesh,
        use_bindings=bindings,
        cbuffer=cbuffer,
    )

    rpath = str(tmpdir) + "/"

    assert mesh.validate(rpath=rpath, scribes=0)
