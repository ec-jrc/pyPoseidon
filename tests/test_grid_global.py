import pyPoseidon.grid as pg
import numpy as np
import pytest
import os
import geopandas as gp
import cartopy.feature as cf


ne=[]
for cr in ['l','i','h']:
    coast = cf.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='{}m'.format({'l':110, 'i':50, 'h':10}[cr]))

    gdf = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])

    w = gdf.explode().reset_index(drop=True)
    ne.append(w)

coast = cf.GSHHSFeature(
    scale='auto',
    levels = [1])

GSHHS = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])


@pytest.mark.slow
@pytest.mark.parametrize('coast', ne[0:1])

def test_answer(tmpdir, coast):

    df = pg.grid(type='tri2d',geometry='global', coastlines=coast, rpath = str(tmpdir)+'/')

    check = np.isnan(df.Dataset.depth.values).sum() == 0

    assert check == True
