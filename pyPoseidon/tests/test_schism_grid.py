import pyPoseidon.grid as pgrid
import pytest
import os


def func(tmpdir,name):

    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', name)
    #read grid file
    grid = pgrid.grid(type='tri2d',grid_file=filename)

    filename_ = str(tmpdir.join('hgrid_.gr3'))
    #output to grid file
    grid.impl.to_file(filename_)

    #read again new grid
    grid_ = pgrid.grid(type='tri2d',grid_file=filename_)

    #cleanup
    os.remove(filename_)

    #compare
    return grid.impl.Dataset.equals(grid_.impl.Dataset)


def test_answer(tmpdir):
    assert func(tmpdir,'hgrid.gr3') == True