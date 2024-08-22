from __future__ import annotations

import numpy as np
import pandas as pd
import sklearn.neighbors
from scipy.spatial.distance import cdist


def closest_node(node, nodes):
    return nodes[cdist([node], nodes).argmin()]


def find_nearest_nodes(
    mesh_nodes: pd.DataFrame,
    points: pd.DataFrame,
    metric: str = "haversine",
    earth_radius = 6371000,
    ):
    """
    Calculate the mesh nodes that are nearest to the specified `points`.

    Both `mesh_nodes` and `points` must be `pandas.DataFrames` that have
    columns named `lon` and `lat` and the coords must be in EPSG:4326.

    Returns the `points` DataFrame after adding these extra columns:

    - `mesh_index` which is the index of the node in the `hgrid.gr3` file
    - `mesh_lon` which is the longitude of the nearest mesh node
    - `mesh_lat` which is the latitude of the nearest mesh node
    - `distance` which is the distance in meters between the point and the nearest mesh node

    Examples:
        >>> mesh_nodes = pd.DataFrame({
        ...     "lon": [0, 10, 20],
        ...     "lat": [0, 5, 0],
        ... })
        >>> points = pd.DataFrame({
        ...     "lon": [1, 11, 21],
        ...     "lat": [1, 4, 1],
        ...     "id": ["a", "b", "c"],
        ... })
        >>> nearest_nodes = find_nearest_nodes(mesh_nodes, points)
        >>> nearest_nodes
           lon  lat id  mesh_index  mesh_lon  mesh_lat       distance
        0    1    1  a           0         0         0  157249.381272
        1   11    4  b           1        10         5  157010.162641
        2   21    1  c           2        20         0  157249.381272
    """
    # The only requirement is that both `mesh_nodes and `points` have `lon/lat` columns
    tree = sklearn.neighbors.BallTree(
        np.radians(mesh_nodes[["lat", "lon"]]),
        metric=metric,
    )
    distances, indices = tree.query(np.radians(points[["lat", "lon"]].values))
    closest_nodes = (
        mesh_nodes
        .rename(columns={"lon": "mesh_lon", "lat": "mesh_lat"})
        .iloc[indices.flatten()]
        .assign(distance=(distances.flatten() * earth_radius))
        .reset_index(names=["mesh_index"])
    )
    return pd.concat((points, closest_nodes), axis="columns")
