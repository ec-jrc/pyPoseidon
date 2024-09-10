from __future__ import annotations

import pandas as pd
import pytest

from pyposeidon.utils.cpoint import find_nearest_nodes
from pyposeidon.utils.cpoint import get_ball_tree

EXPECTED_COLUMNS = ["lon", "lat", "id", "mesh_index", "mesh_lon", "mesh_lat", "distance"]


@pytest.fixture(scope="session")
def mesh_nodes():
    return pd.DataFrame({
        "lon": [0, 10, 20],
        "lat": [0, 5, 0],
    })


@pytest.fixture(scope="session")
def points():
    return pd.DataFrame({
        "lon": [1, 11, 21, 2],
        "lat": [1, 4, 1, 2],
        "id": ["a", "b", "c", "d"],
    })


@pytest.fixture(scope="session")
def ball_tree(mesh_nodes):
    return get_ball_tree(mesh_nodes)


def test_find_nearest_nodes(mesh_nodes, points):
    nearest_nodes = find_nearest_nodes(mesh_nodes, points)
    assert isinstance(nearest_nodes, pd.DataFrame)
    assert len(nearest_nodes) == len(points)
    assert nearest_nodes.columns.tolist() == EXPECTED_COLUMNS
    assert nearest_nodes.mesh_index.tolist() == [0, 1, 2, 0]
    assert nearest_nodes.distance.min() > 150_000
    assert nearest_nodes.distance.max() < 320_000


@pytest.mark.parametrize("k", [pytest.param(2, id='2 points'), pytest.param(3, id='3 points')])
def test_find_nearest_nodes_multiple_points_and_pass_tree_as_argument(mesh_nodes, points, k, ball_tree):
    nearest_nodes = find_nearest_nodes(mesh_nodes, points, k=k, tree=ball_tree)
    assert isinstance(nearest_nodes, pd.DataFrame)
    assert len(nearest_nodes) == len(points) * k
    assert nearest_nodes.columns.tolist() == EXPECTED_COLUMNS
