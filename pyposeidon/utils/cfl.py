from __future__ import annotations

import io
import os

import numpy as np
import numpy.typing as npt
import pandas as pd
import pymap3d


def parse_hgrid_nodes(path: os.PathLike[str] | str) -> pd.DataFrame:
    with open(path, "rb") as fd:
        _ = fd.readline()
        _, no_points = map(int, fd.readline().strip().split(b" "))
        content = io.BytesIO(b"".join(next(fd) for _ in range(no_points)))
        nodes = pd.read_csv(
            content,
            engine="pyarrow",
            sep="\t",
            header=None,
            names=["lon", "lat", "depth"],
            index_col=0,
        )
    nodes = nodes.reset_index(drop=True)
    return nodes


def parse_hgrid_elements3(path: os.PathLike[str] | str) -> pd.DataFrame:
    with open(path, "rb") as fd:
        _ = fd.readline()
        no_elements, no_points = map(int, fd.readline().strip().split(b" "))
        for _ in range(no_points):
            next(fd)
        content = io.BytesIO(b"".join(next(fd) for _ in range(no_elements)))
        elements = pd.read_csv(
            content,
            engine="pyarrow",
            sep="\t",
            header=None,
            names=["no_nodes", "n1", "n2", "n3"],
            index_col=0,
        )
    elements = elements.assign(
        n1=elements.n1 - 1,
        n2=elements.n2 - 1,
        n3=elements.n3 - 1,
    ).reset_index(drop=True)
    return elements


def get_skews_and_base_cfls(lons, lats, depths) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
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
    base_cfls = np.sqrt(9.81 * np.maximum(0.1, depths.mean(axis=0))) / rhos / 2
    return skews, base_cfls


def get_skews_and_base_cfls_from_path(
    path: os.PathLike[str] | str,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    nodes = parse_hgrid_nodes(path)
    elements = parse_hgrid_elements3(path)
    tri = elements[["n1", "n2", "n3"]].values
    lons = nodes.lon.values[tri].T
    lats = nodes.lat.values[tri].T
    depths = nodes.depth.values[tri].T
    skews, base_cfls = get_skews_and_base_cfls(lons=lons, lats=lats, depths=depths)
    return skews, base_cfls
