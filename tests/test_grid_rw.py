import pyPoseidon.grid as pg
import numpy as np
import pytest
import os

from . import DATA_DIR

GRID_FILE = DATA_DIR / "hgrid.gr3"

def test_answer(tmpdir):
    
    rpath = str(tmpdir)+'/'
    
    
    df = pg.grid(type='tri2d',grid_file=GRID_FILE)
    
    df.to_file(rpath + 'test.gr3')
    
    dh = pg.grid(type='tri2d',grid_file=rpath+'test.gr3')
    
    assert df.Dataset.equals(dh.Dataset)