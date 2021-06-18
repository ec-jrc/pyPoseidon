import pyposeidon.grid as pg
import numpy as np
import pytest
import os
import geopandas as gp
import cartopy.feature as cf

from . import DATA_DIR

DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

COAST_FILE = (DATA_DIR / "ocean.zip").as_posix()


@pytest.mark.slow
@pytest.mark.parametrize("ggor", ["jigsaw", "gmsh"])
@pytest.mark.parametrize("bgmesh", [None, DEM_FILE])
@pytest.mark.parametrize("bindings", [True, False])
def test_answer(tmpdir, ggor, bgmesh, bindings):

    df = pg.grid(
        type="tri2d",
        geometry="global",
        coastlines=COAST_FILE,
        rpath=str(tmpdir) + "/",
        grid_generator=ggor,
        dem_source=bgmesh,
        use_bindings=bindings,
    )

    check = np.isnan(df.Dataset.depth.values).sum() == 0

    assert check == True
