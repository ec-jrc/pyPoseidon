from __future__ import annotations

import pandas as pd
from pyposeidon.utils.cpoint import find_nearest_nodes


def test_find_nearest_nodes():
    mesh_nodes = pd.DataFrame({
        "lon": [0, 10, 20],
        "lat": [0, 5, 0],
    })
    points = pd.DataFrame({
        "lon": [1, 11, 21],
        "lat": [1, 4, 1],
        "id": ["a", "b", "c"],
    })
    nearest_nodes = find_nearest_nodes(mesh_nodes, points)
    assert isinstance(nearest_nodes, pd.DataFrame)
    assert len(nearest_nodes) == 3
    assert nearest_nodes.columns.tolist() == ["lon", "lat", "id", "mesh_index", "mesh_lon", "mesh_lat", "distance"]
    assert nearest_nodes.mesh_index.tolist() == [0, 1, 2]
    assert nearest_nodes.distance.min() > 150_000
    assert nearest_nodes.distance.max() < 160_000
