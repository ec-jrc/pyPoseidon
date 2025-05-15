from __future__ import annotations

import collections
import io
import itertools
import os
import typing as T

import numpy as np
import numpy.typing as npt
import pymap3d


def _readline(fd: bytes) -> bytes:
    return fd.readline().split(b"=")[0].split(b"!")[0].strip()


def parse_hgrid(
    path: os.PathLike[str] | str,
    include_boundaries: bool = False,
    sep: str | None = None,
) -> dict[str, T.Any]:
    """
    Parse an hgrid.gr3 file.

    The function is also able to handle fort.14 files, too, (i.e. ADCIRC)
    but the boundary parsing is not keeping all the available information.
    """
    rvalue: dict[str, T.Any] = {}
    with open(path, "rb") as fd:
        _ = fd.readline()  # skip line
        no_elements, no_points = map(int, fd.readline().strip().split(b"!")[0].split())
        nodes_buffer = io.BytesIO(b"\n".join(itertools.islice(fd, 0, no_points)))
        nodes = np.loadtxt(nodes_buffer, delimiter=sep, usecols=(1, 2, 3))
        elements_buffer = io.BytesIO(b"\n".join(itertools.islice(fd, 0, no_elements)))
        elements = np.loadtxt(elements_buffer, delimiter=sep, usecols=(2, 3, 4), dtype=int)
        elements -= 1  # 0-based index for the nodes
        rvalue["nodes"] = nodes
        rvalue["elements"] = elements
        # boundaries
        if include_boundaries:
            boundaries = collections.defaultdict(list)
            no_open_boundaries = int(_readline(fd))
            total_open_boundary_nodes = int(_readline(fd))
            for i in range(no_open_boundaries):
                no_nodes_in_boundary = int(_readline(fd))
                boundary_nodes = np.loadtxt(fd, delimiter=sep, usecols=(0,), dtype=int)
                boundaries["open"].append(boundary_nodes - 1)  # 0-based index
            # closed boundaries
            no_closed_boundaries = int(_readline(fd))
            total_closed_boundary_nodes = int(_readline(fd))
            for _ in range(no_closed_boundaries):
                # Sometimes it seems that the closed boundaries don't have a "type indicator"
                # For example: Test_COSINE_SFBay/hgrid.gr3
                # In this cases we assume that boundary type is 0 (i.e. land in schism)
                # XXX Maybe check the source code?
                parsed = _readline(fd).split(b" ")
                if len(parsed) == 1:
                    no_nodes_in_boundary = int(parsed[0])
                    boundary_type = 0
                else:
                    no_nodes_in_boundary, boundary_type = map(int, (p for p in parsed if p))
                boundary_nodes = np.genfromtxt(fd, delimiter=sep, usecols=(0,), max_rows=no_nodes_in_boundary, dtype=int)
                boundary_nodes -= 1  # 0-based-index
                boundaries[boundary_type].append(boundary_nodes)
            rvalue["boundaries"] = boundaries
    return rvalue


def get_skews_and_base_cfls(
    lons,
    lats,
    depths,
    g: float = 9.81,
    minimum_depth: float = 0.1,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    # The shape of each one of the input arrays needs to be (3, <no_triangles>)
    # ell = pymap3d.Ellipsoid.from_name("wgs84")
    ell = pymap3d.Ellipsoid(6378206.4, 6378206.4, "schism", "schism")
    local_x, local_y, _ = pymap3d.geodetic2enu(lats, lons, depths, lats[0], lons[0], depths[0], ell=ell)
    areas = (local_x[1] * local_y[2] - local_x[2] * local_y[1]) * 0.5
    rhos = np.sqrt(areas / np.pi)
    max_sides = np.maximum(
        np.sqrt(local_x[1] ** 2 + local_y[1] ** 2),
        np.sqrt(local_x[2] ** 2 + local_y[2] ** 2),
        np.sqrt((local_x[2] - local_x[1]) ** 2 + (local_y[2] - local_y[1]) ** 2),
    )
    skews = max_sides / rhos
    base_cfls = np.sqrt(g * np.maximum(minimum_depth, depths.mean(axis=0))) / rhos / 2
    return skews, base_cfls

def get_skews_and_base_cfls_from_path(
    path: os.PathLike[str] | str,
    g: float = 9.81,
    minimum_depth: float = 0.1,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    parsed = parse_hgrid(path)
    tri = parsed["elements"]
    nodes = parsed["nodes"]
    lons = nodes[:, 0][tri].T
    lats = nodes[:, 1][tri].T
    depths = nodes[:, 2][tri].T
    skews, base_cfls = get_skews_and_base_cfls(lons=lons, lats=lats, depths=depths, g=g, minimum_depth=minimum_depth)
    return skews, base_cfls
