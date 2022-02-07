import pyposeidon.mesh as pmesh
import pytest
import os

from . import DATA_DIR


def func(tmpdir, name):

    filename = str(DATA_DIR / name)
    # read mesh file
    mesh = pmesh.set(type="tri2d", mesh_file=filename)

    filename_ = str(tmpdir.join("hgrid_.gr3"))
    # output to mesh file
    mesh.to_file(filename_)

    # read again new mesh
    mesh_ = pmesh.set(type="tri2d", mesh_file=filename_)

    # cleanup
    os.remove(filename_)

    # compare
    return mesh.Dataset.equals(mesh_.Dataset)


def test_answer(tmpdir):
    assert func(tmpdir, "hgrid.gr3") == True
