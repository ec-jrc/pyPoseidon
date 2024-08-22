import pyposeidon.mesh as pmesh
import numpy as np
import pytest
import os

from . import DATA_DIR

COAST_FILE = (DATA_DIR / "ocean.parquet").as_posix()

WINDOWS = pytest.mark.parametrize(
    "window",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 165.0, "lon_max": 185.0, "lat_min": -50, "lat_max": -30},
        {"lon_min": -185.0, "lon_max": -165.0, "lat_min": -50, "lat_max": -30},
        {"lon_min": -25.0, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 68.0},
    ],
)


@WINDOWS
@pytest.mark.parametrize("mesh_generator", ["jigsaw", "gmsh"])
def test_answer(tmpdir, window, mesh_generator):
    mesh = pmesh.set(
        type="tri2d",
        geometry=window,
        coastlines=COAST_FILE,
        rpath=str(tmpdir) + "/",
        mesh_generator=mesh_generator,
    )

    check = np.isnan(mesh.Dataset.depth.values).sum() == 0
    assert check == True
