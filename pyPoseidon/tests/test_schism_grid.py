import pyPoseidon.grid as pgrid
import pytest
import os

from . import DATA_DIR


def func(tmpdir,name):

    filename = DATA_DIR / name
    #read grid file
    grid = pgrid.grid(type='tri2d',grid_file=filename)

    filename_ = str(tmpdir.join('hgrid_.gr3'))
    #output to grid file
    grid.to_file(filename_)

    #read again new grid
    grid_ = pgrid.grid(type='tri2d',grid_file=filename_)

    #cleanup
    os.remove(filename_)

    #compare
    return grid.Dataset.equals(grid_.Dataset)


def test_answer(tmpdir):
    assert func(tmpdir,'hgrid.gr3') == True
