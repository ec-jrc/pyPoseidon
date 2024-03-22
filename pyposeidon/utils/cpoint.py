from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
import numpy as np


def closest_node(node, nodes):
    return nodes[cdist([node], nodes).argmin()]


def closest_n_points(nodes, N, meshXY, dist_max=np.inf):
    """
    Find the indices of the N closest points in a set of points.

    Parameters
    ----------
    nodes : np.ndarray
        The set of points to search.
    N : int
        The number of points to find.
    meshXY : np.ndarray
        The grid of points used for the KDTree.
    dist_max : float, optional (default np.inf)
        the max distance to reach to interpolate from the nearest point.
        above that limit, no interpolation is done. the point is removed.

    Returns
    -------
    np.ndarray
        The indices of the N closest points.
    np.boolarray
        The mask of the retained stations (with dist_max filter)
    """

    mytree = cKDTree(meshXY)
    dist, indice = mytree.query(nodes, range(1, N + 1))
    indice[dist > dist_max] = -1
    mask = indice != -1
    return indice[mask].T, mask.T[0]
