import sys

import numpy as np
import pytest

import pyposeidon.mesh as pmesh
from . import DATA_DIR

DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

COAST_FILE = (DATA_DIR / "ocean.parquet").as_posix()


@pytest.mark.parametrize("mesh_generator,use_bindings", [("jigsaw", None), ("gmsh", True), ("gmsh", False)])
@pytest.mark.parametrize("dem_source", [None, DEM_FILE])
@pytest.mark.parametrize("cbuffer", [None, 0.01])
def test_io(pytestconfig, tmpdir, mesh_generator, use_bindings, dem_source, cbuffer):
    # Skip the test unless --runslow has been passed
    if dem_source is not None:
        if not pytestconfig.getoption("--runslow"):
            pytest.skip("slow test")

    if mesh_generator == "jigsaw" and cbuffer is not None and sys.platform != "darwin":
        pytest.xfail("jigsaw + buffer is failing on linux: https://github.com/ec-jrc/pyPoseidon/issues/84")

    mesh = pmesh.set(
        type="tri2d",
        geometry="global",
        coastlines=COAST_FILE,
        rpath=str(tmpdir) + "/",
        mesh_generator=mesh_generator,
        dem_source=dem_source,
        use_bindings=use_bindings,
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
@pytest.mark.parametrize("mesh_generator,use_bindings", [("jigsaw", None), ("gmsh", True), ("gmsh", False)])
@pytest.mark.parametrize("dem_source", [None, DEM_FILE])
@pytest.mark.parametrize("cbuffer", [None, 0.01])
def test_validate(pytestconfig, tmpdir, mesh_generator, use_bindings, dem_source, cbuffer):
    if dem_source is not None:
        if not pytestconfig.getoption("--runslow"):
            pytest.skip("slow test")

    if mesh_generator == "jigsaw" and cbuffer is not None and sys.platform != "darwin":
        pytest.xfail("jigsaw + buffer is failing on linux: https://github.com/ec-jrc/pyPoseidon/issues/84")

    mesh = pmesh.set(
        type="tri2d",
        geometry="global",
        coastlines=COAST_FILE,
        rpath=str(tmpdir) + "/",
        mesh_generator=mesh_generator,
        dem_source=dem_source,
        use_bindings=use_bindings,
        cbuffer=cbuffer,
    )

    rpath = str(tmpdir) + "/"

    assert mesh.validate(rpath=rpath, scribes=0)
