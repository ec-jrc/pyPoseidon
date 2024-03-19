from __future__ import annotations

import typing as T

import numpy as np
import pandas as pd
import sklearn.neighbors
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from scipy.spatial import Delaunay
import numpy as np


def closest_node(node, nodes):
    return nodes[cdist([node], nodes).argmin()]


def get_ball_tree(
    mesh_nodes: pd.DataFrame,
    metric: str = "haversine",
) -> sklearn.neighbors.BallTree:
    """
    Return a `BallTree` constructed from the provided `mesh_nodes`.

    `mesh_nodes` must be a `pandas.DataFrames` with columns named
    `lon` and `lat` and the coords must be in EPSG:4326.
    """
    tree = sklearn.neighbors.BallTree(
        np.radians(mesh_nodes[["lat", "lon"]]),
        metric=metric,
    )
    return tree


def find_nearest_nodes(
    mesh_nodes: pd.DataFrame,
    points: pd.DataFrame,
    k: int = 1,
    metric: str = "haversine",
    earth_radius=6371000,
    tree: sklearn.neighbors.BallTree | None = None,
    **kwargs: T.Any,
):
    """
    Calculate the `k` mesh nodes that are nearest to the specified `points`.

    Both `mesh_nodes` and `points` must be `pandas.DataFrames` that have
    columns named `lon` and `lat` and the coords must be in EPSG:4326.

    As a speed optimization, the `tree` can be pre-constructed with ``get_ball_tree()``
    (and/or serialized to disk using [skops](skops.io) or `pickle`)
    and passed using the `tree` argument.

    Returns the `points` DataFrame after adding these extra columns:

    - `mesh_index` which is the index of the mesh node
    - `mesh_lon` which is the longitude of the nearest mesh node
    - `mesh_lat` which is the latitude of the nearest mesh node
    - `distance` which is the distance in meters between the point in question
       and the nearest mesh node

    Examples:
        ``` python
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
        >>> nearest_nodes = find_nearest_nodes(mesh_nodes, points, k=2)
        >>> nearest_nodes
           lon  lat id  mesh_index  mesh_lon  mesh_lat       distance
        0    1    1  a           0         0         0  1.572494e+05
        1    1    1  a           1        10         5  1.093700e+06
        2   11    4  b           1        10         5  1.570102e+05
        3   11    4  b           2        20         0  1.094398e+06
        4   21    1  c           2        20         0  1.572494e+05
        5   21    1  c           1        10         5  1.299688e+06
        ```
    """
    # Resolve tree
    if tree is None:
        tree = get_ball_tree(mesh_nodes=mesh_nodes, metric=metric)
    distances, indices = tree.query(
        X=np.radians(points[["lat", "lon"]].values),
        k=k,
        return_distance=True,
        **kwargs,
    )
    closest_nodes = (
        mesh_nodes.rename(columns={"lon": "mesh_lon", "lat": "mesh_lat"})
        .iloc[indices.flatten()]
        .assign(distance=(distances.flatten() * earth_radius))
        .reset_index(names=["mesh_index"])
    )
    return pd.concat((points.loc[points.index.repeat(k)].reset_index(drop=True), closest_nodes), axis="columns")


def get_weights(in_xy, out_xy, d=2):
    """
    Triangulate an output mesh based on a set of input points and return the
    corresponding barycentric weights.

    Parameters
    ----------
    @param in_xy (np.ndarray): Nx2 array of input point coordinates.
    @param out_xy (np.ndarray): Nx2 array of output point coordinates.
    @param d (int): Dimension of the input points, by default 2.

    @returns (4-uple): (vert, wgts, out_idx_out, in_idx_out), with:
        vert (np.ndarray): Nx3 array of triangle vertices.
        wgts (np.ndarray): Nx3 array of barycentric weights.
        out_idx_out (np.ndarray): Boolean array indicating which points
            are outside the input mesh.
        in_idx_out (np.ndarray): Index array of the nearest input points
            for points outside the input mesh.

    """
    t = Delaunay(in_xy)  # triangulate output mesh
    s = t.find_simplex(out_xy)
    vert = np.take(t.simplices, np.maximum(s, 0), axis=0)  # Use max to avoid negative indices
    t_ = np.take(t.transform, np.maximum(s, 0), axis=0)
    delta = out_xy - t_[:, d]
    bary = np.einsum("njk,nk->nj", t_[:, :d, :], delta)
    wgts = np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))
    # Points outside the out_xy
    out_idx_out = s < 0
    if np.any(out_idx_out):
        # For points outside, find nearest neighbors
        tree = cKDTree(in_xy)
        _, in_idx_out = tree.query(out_xy[out_idx_out])
    else:
        in_idx_out = None
    return vert, wgts, out_idx_out, in_idx_out


def interp(values, vtx, wts, out_idx_out, in_idx_out):
    """
    Interpolate values using barycentric coordinates and weights

    Parameters
    ----------
    @param values (np.ndarray): Array of values to be interpolated.
    @param vtx (np.ndarray): Array of vertex indices.
    @param wts (np.ndarray): Array of barycentric weights.
    @param out_idx_out (np.ndarray): Boolean array indicating which points
        are outside the input mesh.
    @param in_idx_out (np.ndarray): Index array of the nearest input points
        for points outside the input mesh.

    @return (np.ndarray): Interpolated values

    """
    res = np.einsum("nj,nj->n", np.take(values, vtx), wts)
    if in_idx_out is not None:
        res[out_idx_out] = values[in_idx_out]
    return res
