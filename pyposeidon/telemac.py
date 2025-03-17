"""
TELEMAC model of pyposeidon. It controls the creation, execution & output  of a complete simulation based on SCHISM.

"""

# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import os
import datetime
import numpy as np
import sys
import json
import pandas as pd
import glob
import xarray as xr
import geopandas as gp
import shapely
import jinja2
from searvey import ioc
import shutil
import typing as T
from scipy.spatial import Delaunay
import dask.array as da

# local modules
import pyposeidon
import pyposeidon.mesh as pmesh
import pyposeidon.meteo as pmeteo
import pyposeidon.dem as pdem
from pyposeidon.paths import DATA_PATH
from pyposeidon.utils.get_value import get_value
from pyposeidon.utils.converter import myconverter
from pyposeidon.utils.cpoint import find_nearest_nodes
from pyposeidon.utils.cpoint import get_weights
from pyposeidon.utils.cpoint import interp
from pyposeidon.utils import data
from pyposeidon.utils.norm import normalize_column_names
from pyposeidon.utils.norm import normalize_varnames
from pyposeidon.utils.obs import serialize_stations
from pyposeidon.utils.post import export_xarray
from pyposeidon.utils.post import remove
from pyposeidon.utils.cast import copy_files

# telemac (needs telemac stack on conda)
from telapy.api.t2d import Telemac2d
from telapy.api.t3d import Telemac3d
from telapy.api.art import Artemis
from telapy.api.wac import Tomawac
from execution.telemac_cas import TelemacCas
from data_manip.formats.selafin import Selafin
from data_manip.extraction.telemac_file import TelemacFile
from utils.progressbar import ProgressBar
from pretel.manip_telfile import spherical2longlat
from pretel.manip_telfile import longlat2spherical
from pretel.extract_contour import extract_boundaries
from pretel.extract_contour import sort_contours
from pretel.extract_contour import signed_area
from collections import defaultdict
from xarray_selafin.xarray_backend import SelafinAccessor
from mpi4py import MPI

import numpy as np

import logging

from . import tools

# Types
ResultTypes = T.Literal["1D", "2D", "3D"]

logger = logging.getLogger(__name__)


TELEMAC_NAME = "telemac"


def divisors(n):
    result = set()
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            result.add(i)
            result.add(n // i)
    return sorted(result)


def calculate_time_step_hourly_multiple(resolution, min_water_depth=0.1, courant_number=50):
    """
    Calculate the maximum allowable time step for a shallow water equation model
    based on the CFL condition, rounded to the closest divisor of an hour.

    :param resolution: Spatial resolution (Δx) in meters
    :param min_water_depth: Minimum water depth (h_min) in meters. Default is 1 meter.
    :param courant_number: Courant number (C), typically less than 1 for explicit schemes. Can be higher for implicit schemes.
    :return: Maximum allowable time step (Δt) in seconds, rounded to the closest divisor of an hour.
    """
    resolution *= 111.389  # degrees to meters conversion factor
    raw_time_step = courant_number * resolution / (9.81 * min_water_depth) ** 0.5
    # Get all divisors of an hour (3600 seconds)
    hour_divisors = divisors(3600)
    # Find the closest divisor of an hour to the raw time step
    closest_divisor = min(hour_divisors, key=lambda x: abs(x - raw_time_step))
    return closest_divisor


# Helper functions
def df_to_gpd(df, x: str = "longitude", y: str = "latitude"):
    gpdf = gp.GeoDataFrame(
        df,
        geometry=gp.points_from_xy(
            df[x],
            df[y],
        ),
        crs="EPSG:4326",
    )
    return gpdf


def xy_to_ll(x, y):
    gpdf = gp.GeoDataFrame(geometry=gp.points_from_xy(x, y), crs="EPSG:3857").to_crs("EPSG:4326")
    return gpdf.geometry.x, gpdf.geometry.y


def is_ccw(tris, meshx, meshy):
    x1, x2, x3 = meshx[tris].T
    y1, y2, y3 = meshy[tris].T
    return (y3 - y1) * (x2 - x1) > (y2 - y1) * (x3 - x1)


def flip(tris):
    return np.column_stack((tris[:, 2], tris[:, 1], tris[:, 0]))


def is_overlapping(tris, meshx):
    PIR = 180
    x1, x2, x3 = meshx[tris].T
    return np.logical_or(abs(x2 - x1) > PIR, abs(x3 - x1) > PIR, abs(x3 - x3) > PIR)


def get_det_mask(tris, meshx, meshy):
    t12 = meshx[tris[:, 1]] - meshx[tris[:, 0]]
    t13 = meshx[tris[:, 2]] - meshx[tris[:, 0]]
    t22 = meshy[tris[:, 1]] - meshy[tris[:, 0]]
    t23 = meshy[tris[:, 2]] - meshy[tris[:, 0]]
    return t12 * t23 - t22 * t13 > 0


def contains_pole(x, y):
    return np.any(y == 90, axis=0)


def flip_triangles(x, y, ikle2, corrections):
    # Assuming slf.meshx, slf.meshy, and slf.ikle2 are NumPy arrays
    ikle2 = np.array(ikle2)
    # Ensure all triangles are CCW
    ccw_mask = is_ccw(ikle2, x, y)
    ikle2[~ccw_mask] = flip(ikle2[~ccw_mask])

    # triangles accross the dateline
    m_ = is_overlapping(ikle2, x)
    ikle2[m_] = flip(ikle2[m_])

    # manual additional corrections
    if corrections is not None:
        for rev in corrections["reverse"]:
            ikle2[rev : rev + 1] = flip(ikle2[rev : rev + 1])

    # Check for negative determinant
    detmask = ~get_det_mask(ikle2, x, y)
    logger.debug(f"reversed {detmask.sum()} triangles")
    logger.debug("[TEMPORARY FIX]: REMOVE THE POLE TRIANGLES")
    return ikle2


def remove_pole_triangles(x, y, ikle2, corrections):
    # Assuming slf.meshx, slf.meshy, and slf.ikle2 are NumPy arrays
    ikle2 = np.array(ikle2)
    pole_mask = contains_pole(x[ikle2].T, y[ikle2].T)
    # manual additional corrections
    if corrections is not None:
        for rem in corrections["remove"]:
            pole_mask[rem] = True

    logger.debug(f"pole triangles: {np.where(pole_mask)[0]}")
    logger.debug("[TEMPORARY FIX]: REMOVE THE POLE TRIANGLES")
    return ikle2[~pole_mask, :]


def write_netcdf(ds, outpath):
    fileOut = os.path.splitext(outpath)[0] + ".nc"
    ds.to_netcdf(fileOut)


def extract_t_elev_2D(
    ds: xr.Dataset,
    x: float,
    y: float,
    var: str = "elev",
    xstr: str = "longitude",
    ystr: str = "latitude",
    max_dist: float = 1000,
):
    lons, lats = ds[xstr].values, ds[ystr].values
    mesh = pd.DataFrame(np.vstack([x, y]).T, columns=["lon", "lat"])
    points = pd.DataFrame(np.vstack([lons, lats]).T, columns=["lon", "lat"])
    df = find_nearest_nodes(mesh, points, 1)
    df = df[df.distance < max_dist]
    indx = df["mesh_index"]
    ds_ = ds.isel(node=indx.values[0])
    out_ = ds_[var].values
    t_ = [pd.Timestamp(ti) for ti in ds_.time.values]
    return pd.Series(out_, index=t_), float(ds_[xstr]), float(ds_[ystr])


# Function to subset ERA data based on the mesh extent
def subset_era_from_mesh(
    era: xr.Dataset,
    mesh: xr.Dataset,
    input360: bool,
    gtype: str,
) -> xr.Dataset:
    """
    Selects a subset of ERA data that overlaps with the provided mesh's geographical extent.

    :param era: The ERA dataset from which to select a subset.
    :param mesh: The mesh dataset defining the geographical extent for the subset selection.
    :param input360: A flag indicating whether the longitude should be adjusted to a 0-360 range.
    :return: A subset of the ERA dataset that falls within the mesh's geographical extent.
    """
    xmin, xmax, ymin, ymax = mesh.x.min(), mesh.x.max(), mesh.y.min(), mesh.y.max()
    if input360:
        xmin, xmax = np.mod(xmin + 360, 360), np.mod(xmax + 360, 360)
        if xmax <= xmin:
            xmin, xmax = 0, 360
    if gtype == "grid":
        era_chunk = era.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))
    else:  # for 01280 grid
        mask_lon = (era.longitude >= xmin) & (era.longitude <= xmax)
        mask_lat = (era.latitude >= ymin) & (era.latitude <= ymax)
        mask = mask_lon & mask_lat
        indices = np.where(mask)[0]
        era_chunk = era.isel(values=indices)
    return era_chunk


# Function to write meteorological data onto a mesh
def write_meteo_on_mesh(
    era_ds: xr.Dataset,
    mesh: xr.Dataset,
    file_out: str,
    n_time_chunk: int,
    n_node_chunk: int,
    input360: bool = True,
    gtype: str = "grid",
    ttype: str = "time",
) -> None:
    """
    Writes meteorological data from an ERA dataset onto a mesh and saves the result as a zarr file.

    :param era_ds: The ERA dataset with the meteorological data.
    :param mesh: The mesh dataset representing the spatial domain.
    :param file_out: The path to the output zarr file.
    :param n_time_chunk: The size of the time chunks for processing.
    :param n_node_chunk: The size of the node chunks for processing.
    :param input360: A flag indicating whether the longitude should be adjusted to a 0-360 range.
    """
    # Create the temporary dummy zarr file
    if os.path.exists(file_out):
        shutil.rmtree(file_out)
    x, y, tri = mesh.x.values, mesh.y.values, mesh.attrs["ikle2"] - 1
    nnodes = len(x)
    ntimes = len(era_ds.time)
    zero = da.zeros((ntimes, nnodes), chunks=(n_time_chunk, n_node_chunk))

    # Define coordinates and data variables for the output dataset
    coords = {
        "time": era_ds.time,
        "node": np.arange(nnodes),
        "lon": ("node", x),
        "lat": ("node", y),
        "triface_nodes": (("face_nodes", "three"), tri),
    }
    data_vars = {}
    for varin in era_ds.data_vars:
        data_vars[varin] = (("time", "node"), zero)
    xr.Dataset(data_vars=data_vars, coords=coords).to_zarr(file_out, compute=False)

    # in the case of "tstps"
    if ttype == "step":
        t0 = pd.Timestamp(era_ds.time.values)
        seconds = era_ds.step.values / 1e9
        era_ds.time = pd.to_datetime(t0 + pd.Timedelta(seconds=seconds))

    # Iterate over nodes in chunks and write data to zarr file
    for ins in range(0, nnodes, n_node_chunk):
        end_node = min(ins + n_node_chunk, nnodes)
        node_chunk = np.arange(ins, end_node)
        mesh_chunk = mesh.isel(node=slice(ins, end_node))
        era_chunk = subset_era_from_mesh(era_ds, mesh_chunk, input360=input360, gtype=gtype)

        # drop needless coords
        main_coords = [c for c in era_chunk.coords]
        main_dims = [d for d in era_chunk.dims]

        # Get weights for interpolation
        if gtype == "grid":
            nx1d = len(era_chunk.longitude)
            ny1d = len(era_chunk.latitude)
            xx = np.tile(era_chunk.longitude, ny1d).reshape(ny1d, nx1d).T.ravel()
            yy = np.tile(era_chunk.latitude, nx1d)
        else:  # for O1280 grids
            main_dims.extend(["longitude", "latitude"])
            xx = era_chunk.longitude.compute()
            yy = era_chunk.latitude.compute()

        drop_coords = set(main_coords) - set(main_dims)
        era_chunk = era_chunk.drop_vars(drop_coords)

        if input360:
            xx[xx > 180] -= 360
        in_xy = np.vstack((xx, yy)).T
        out_xy = np.vstack((mesh_chunk.x, mesh_chunk.y)).T

        logger.debug("computing weigths for interpolation..")
        vert, wgts, u_x, g_x = get_weights(in_xy, out_xy)

        # Interpolate and write data for each variable and time chunk
        for var_name in era_chunk.data_vars:
            for it_chunk in range(0, ntimes, n_time_chunk):
                t_end = min(it_chunk + n_time_chunk, ntimes)
                time_chunk = era_chunk.time[it_chunk:t_end]
                data_chunk = da.zeros((len(time_chunk), len(node_chunk)), chunks=(n_time_chunk, n_node_chunk))
                for it, t_ in enumerate(time_chunk):
                    tmp = np.ravel(np.transpose(era_chunk.isel(time=it_chunk + it)[var_name].values))
                    data_chunk[it, :] = interp(tmp, vert, wgts, u_x, g_x)
                coords = {"time": time_chunk, "node": node_chunk}
                ds = xr.Dataset({var_name: (("time", "node"), data_chunk)}, coords=coords)
                region = {"time": slice(it_chunk, t_end), "node": slice(ins, end_node)}
                ds.to_zarr(file_out, mode="a", region=region)


def write_meteo_selafin(outpath, input_atm_zarr):
    xatm = xr.open_dataset(input_atm_zarr, engine="zarr")
    t0 = pd.Timestamp(xatm.time.values[0])
    # Define a mapping from the original variable names to the new ones
    var_map = {
        "u10": ("WINDX", "M/S"),
        "v10": ("WINDY", "M/S"),
        "msl": ("PATM", "PASCAL"),
        "tmp": ("TAIR", "DEGREES C"),
    }
    var_attrs = {}
    for var in xatm.data_vars:
        if var in var_map:
            # Attributes for the variable
            var_attrs[var] = (var_map[var][0], var_map[var][1])
    # Add global attributes after concatenation
    xatm.attrs["date_start"] = [t0.year, t0.month, t0.day, t0.hour, t0.minute, t0.second]
    xatm.attrs["ikle2"] = xatm.triface_nodes.values + 1
    xatm.attrs["variables"] = {var: attrs for var, attrs in var_attrs.items()}
    xatm = xatm.rename({"lon": "x", "lat": "y"})
    xatm = xatm.drop_vars(["triface_nodes"])
    xatm.selafin.write(outpath)
    remove(input_atm_zarr)


def get_boundary_settings(boundary_type, glo_node, bnd_node, module):
    if module in ["telemac2d", "telemac3d", "tomawac"]:
        settings = {
            "lihbor": 5 if boundary_type == "open" else 2,
            "liubor": 6 if boundary_type == "open" else 2,
            "livbor": 6 if boundary_type == "open" else 2,
            "hbor": 0.0,
            "ubor": 0.0,
            "vbor": 0.0,
            "aubor": 0.0,
            "litbor": 5 if boundary_type == "open" else 2,
            "tbor": 0.0,
            "atbor": 0.0,
            "btbor": 0.0,
            "nbor": glo_node + 1,
            "k": bnd_node + 1,
        }
    return settings


def format_value(value, width, precision=3, is_float=False):
    if is_float:
        return f"{value:{5}.{precision}f}"
    else:
        return f"{value:{2}}"


def get_first_point(ds, bnd_points):
    """
    Determine the southwest points among points given trough their indexes

    @param ds (xr.Dataset) dataset files
    @param bnd_points (numpy array of shape (number of points)) contains indexes
    of points
    @return first_bnd_pt_index (integer) index of the southwest point
    """

    x_plus_y = ds.node_x[bnd_points] + ds.node_y[bnd_points]

    southwest_bnd_pts_index = bnd_points[
        np.where(x_plus_y == x_plus_y.min())[0]
    ]

    if southwest_bnd_pts_index.shape[0] == 1:
        first_bnd_pt_index = southwest_bnd_pts_index[0]
    else:
        first_bnd_pt_index = southwest_bnd_pts_index[
            np.where(
                ds.node_x[southwest_bnd_pts_index]
                == ds.node_x[southwest_bnd_pts_index].min()
            )[0][0]
        ]

    return first_bnd_pt_index


def extract_contour(ds: xr.Dataset):
    """
    Generic function for extraction of contour from a mesh
    @param ds (xr.Dataset) dataset files
    @returns (list) List of polygons
    """
    vertices = np.vstack((ds.node_x, ds.node_y)).T

    # Extract boundary edges and nodes
    boundary_edges = extract_boundaries(ds.face_nodes.values)
    boundary_nodes = set(np.unique(np.hstack(boundary_edges)))

    # Create a list of the neighbours of each node for quick search
    node_neighbours = defaultdict(set)
    for edge in boundary_edges:
        node_neighbours[edge[0]].add(edge[1])
        node_neighbours[edge[1]].add(edge[0])

    # Group boundary nodes coordinates into contours
    contours = []
    contours_idx = []
    while len(boundary_nodes) > 0:
        next_vertex = get_first_point(ds, np.array(list(boundary_nodes)))
        boundary_nodes.remove(next_vertex)
        contour_idx = [next_vertex]
        contour = [vertices[next_vertex].tolist()]
        while True:
            neighbours = node_neighbours[next_vertex].intersection(
                boundary_nodes
            )
            if len(neighbours) == 0:
                break
            next_vertex = neighbours.pop()
            boundary_nodes.remove(next_vertex)
            contour.append(vertices[next_vertex].tolist())
            contour_idx.append(next_vertex)

        # Ensure the contour is closed and append it to the list of contours
        contour.append(contour[0])
        contour_idx.append(contour_idx[0])
        contours.append(contour)
        contours_idx.append(contour_idx)

    # Build the list of domains while ensuring proper contours orientation
    sorted_contours = sort_contours(contours)
    domains = []
    domains_idx = []
    for outer, inners in sorted_contours.items():
        contour = contours[outer]
        contour_idx = contours_idx[outer]
        area = signed_area(contour)
        if area < 0:
            contour = contour[::-1]
            contour_idx = contour_idx[::-1]
        domains.append(contour)
        domains_idx.append(contour_idx)
        if len(inners) > 0:
            for i in inners:
                contour = contours[i]
                contour_idx = contours_idx[i]
                area = signed_area(contour)
                if area > 0:
                    contour = contour[::-1]
                    contour_idx = contour_idx[::-1]
                domains.append(contour)
                domains_idx.append(contour_idx)
        

    return domains, domains_idx

def export_cli(ds: xr.Dataset, tel_path: str, outCli: str, tel_module: str = "telemac2d"):
    """
    (This function is a modification of the existing extract_contour() function
    in scripts/python3/pretel/extract_contour.py of the TELEMAC scripts)

    Generic function for extraction of contour from a mesh (with our without
    boundary file)

    @param ds (xr.Dataset) xarray Dataset of the mesh file (used to extract the boundary types)
    @param geo (str) path to the Telemac geo file
    @param outCli (str) Path to the output contour file

    @returns (list) List of polygons
    """
    # normalise
    ds = ds.rename(
        {
            "SCHISM_hgrid_node_x": "node_x",
            "SCHISM_hgrid_node_y": "node_y",
            "SCHISM_hgrid_face_nodes": "face_nodes",
        }
    )
    #
    node_to_type = dict(zip(ds.bnodes.values, ds.type.values))
    domains_bnd = []
    lines = []
    bnd_node = 0
    contours, contours_idx = extract_contour(ds)
    for bnd in contours_idx:
        poly_bnd = []
        coord_bnd = []
        for i, glo_node in enumerate(bnd[:-1]):  # not taking the last node (not repeating)
            x, y = ds.node_x[glo_node], ds.node_y[glo_node]
            coord_bnd.append((x, y))
            # Determine boundary type for the current node
            boundary_type = node_to_type.get(glo_node, "Unknown")
            if boundary_type == "open":
                # Get previous and next node indices in a circular manner
                prev_node = bnd[i - 1] if i > 0 else bnd[-2]  # -2 to skip the repeated last node
                next_node = bnd[i + 1] if i < len(bnd) - 2 else bnd[0]  # Wrapping around to the first node
                # Get boundary types for previous and next nodes
                prev_boundary_type = node_to_type.get(prev_node, "Unknown")
                next_boundary_type = node_to_type.get(next_node, "Unknown")
                # If both adjacent nodes are not 'open', then bnd is closed
                if prev_boundary_type != "open" and next_boundary_type != "open":
                    boundary_type = "Unknown"
            boundary_settings = get_boundary_settings(boundary_type, glo_node, bnd_node, tel_module)

            keys_order = [
                "lihbor",
                "liubor",
                "livbor",
                "hbor",
                "ubor",
                "vbor",
                "aubor",
                "litbor",
                "tbor",
                "atbor",
                "btbor",
                "nbor",
                "k",
            ]
            line = " ".join(str(boundary_settings[key]) for key in keys_order)
            lines.append(line)
            bnd_node += 1
        poly_bnd.append((coord_bnd, bnd))
    domains_bnd.append(poly_bnd)

    # Writing to file
    with open(outCli, "w") as f:
        for line in lines:
            formatted_line = " ".join(
                format_value(value, 3, is_float=isinstance(value, float)) for value in line.split()
            )
            f.write(f"{formatted_line}\n")

    return domains_bnd


def write_cas(outpath, tel_module, params):
    # Load the Jinja2 template
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.join(DATA_PATH)))

    # Define the custom filter for tidal bounds
    def repeat_value(value, count):
        return ";".join([str(value)] * count)

    env.filters["repeat_value"] = repeat_value
    template = env.get_template(tel_module + ".cas")
    ouput = template.render(params)
    with open(os.path.join(outpath, tel_module + ".cas"), "w") as f:
        f.write(ouput)


class Telemac:
    def __init__(self, **kwargs):
        """
        Create a TELEMAC model:
        by default : TELEMAC2D is used.

        Args:
            rfolder str: the TELEMAC module to use. Can be one of the following: telemac3d, telemac2d, artemis, tomawac


        """
        tel_module = kwargs.get("module", "telemac2d")
        if tel_module not in ["telemac3d", "telemac2d", "artemis", "tomawac"]:
            raise ValueError("module must be one of the following: telemac3d, telemac2d, artemis, tomawac")

        rfolder = kwargs.get("rfolder", None)
        if rfolder:
            self.read_folder(**kwargs)

        self.geometry = kwargs.get("geometry", None)

        if self.geometry:
            if isinstance(self.geometry, dict):
                self.lon_min = self.geometry["lon_min"]
                self.lon_max = self.geometry["lon_max"]
                self.lat_min = self.geometry["lat_min"]
                self.lat_max = self.geometry["lat_max"]
            elif self.geometry == "global":
                logger.warning("geometry is 'global'")
            elif isinstance(self.geometry, str):
                try:
                    geo = gp.GeoDataFrame.from_file(self.geometry)
                except:
                    logger.error("geometry argument not a valid geopandas file")
                    sys.exit(1)

                (
                    self.lon_min,
                    self.lat_min,
                    self.lon_max,
                    self.lat_max,
                ) = geo.total_bounds
            else:
                logger.warning("no geometry given")

        # coastlines
        coastlines = kwargs.get("coastlines", None)

        if coastlines is not None:
            try:
                coast = gp.GeoDataFrame.from_file(coastlines)
            except:
                coast = gp.GeoDataFrame(coastlines)

            self.coastlines = coast

        if not hasattr(self, "start_date"):
            start_date = kwargs.get("start_date", None)
            self.start_date = pd.to_datetime(start_date)

        if not hasattr(self, "end_date"):
            if "time_frame" in kwargs:
                time_frame = kwargs.get("time_frame", None)
                self.end_date = self.start_date + pd.to_timedelta(time_frame)
            elif "end_date" in kwargs:
                end_date = kwargs.get("end_date", None)
                self.end_date = pd.to_datetime(end_date)
                self.time_frame = self.end_date - self.start_date

        # if it is the first computation, set origin time (important for restart)
        if not hasattr(self, "rdate"):
            self.rdate = get_value(self, kwargs, "rdate", self.start_date)

        if not hasattr(self, "end_date"):
            # ---------------------------------------------------------------------
            logger.warning("model not set properly, No end_date\n")
            # ---------------------------------------------------------------------

        # dt is in parameters
        params = kwargs.get("parameters", None)
        if params is not None:
            self.tstep = params.get("dt", None)
        else:
            self.tstep = None

        self.module = tel_module
        self.tide = kwargs.get("tide", False)
        tpxo_ = kwargs.get("tpxo_source", None)
        if tpxo_ is not None:
            self.tide = True

        self.atm = kwargs.get("atm", True)
        self.monitor = kwargs.get("monitor", False)

        self.solver_name = TELEMAC_NAME
        # restart
        restart = get_value(self, kwargs, "hotstart", None)
        if restart:
            self.restart = pd.to_datetime(restart)
        else:
            self.restart = None

        # specific to meteo grib files
        self.gtype = get_value(self, kwargs, "meteo_gtype", "grid")
        self.ttype = get_value(self, kwargs, "meteo_ttype", "time")
        self.ncsize = get_value(self, kwargs, "ncsize", 1)
        # convert -180/180 to 0-360
        self.input360 = get_value(self, kwargs, "meteo_input360", False)
        self.meteo = get_value(self, kwargs, "meteo_source", None)

        # custom user fortran
        self.fortran = get_value(self, kwargs, "fortran", None)

        for attr, value in kwargs.items():
            if not hasattr(self, attr):
                setattr(self, attr, value)

        setattr(self, "misc", {})

    # ============================================================================================
    # CONFIG
    # ============================================================================================
    def config(self, config_file=None, output=False, **kwargs):
        dic = get_value(self, kwargs, "parameters", None)
        #        param_file = get_value(self,kwargs,'config_file',None)
        path = get_value(self, kwargs, "rpath", "./telemac/")

        if config_file:
            # ---------------------------------------------------------------------
            logger.info("reading parameter file {}\n".format(config_file))
            # ---------------------------------------------------------------------
        else:
            # ---------------------------------------------------------------------
            logger.info("using default " + self.module + " json config file ...\n")
            # ---------------------------------------------------------------------

            config_file = os.path.join(DATA_PATH, self.module + ".json")

        # Load the parameters from the JSON file
        with open(config_file, "r") as json_file:
            config = json.load(json_file)
            params = config["params"]

        res_min = get_value(self, kwargs, "resolution_min", 0.5)
        # update key values
        if self.start_date is not None:
            if "telemac" in self.module:
                params["datestart"] = self.start_date.strftime("%Y;%m;%d")
                params["timestart"] = self.start_date.strftime("%H;%M;%S")
            elif "tomawac" in self.module:
                params["datestart"] = self.start_date.strftime("%Y%m%d%H%M")
            else:
                logger.info(f"self.start_date not set for {self.module} yet\n")

            if hasattr(self, "time_frame"):
                duration = pd.to_timedelta(self.time_frame).total_seconds()
            else:
                self.time_frame = self.end_date - self.start_date
                duration = self.time_frame.total_seconds()

            params["duration"] = duration

            # export grid data every hour
            if self.tstep:
                tstep = self.tstep
            else:
                tstep = calculate_time_step_hourly_multiple(res_min)
            params["tstep"] = tstep
            params["nb_tsteps"] = int(duration / tstep)
            params["tstep_graph"] = int(3600 / tstep)
            params["tstep_list"] = int(3600 / tstep)
            params["ncsize"] = self.ncsize

        # tide
        if self.tide:
            params["tide"] = True
            params["initial_conditions"] = "TPXO SATELLITE ALTIMETRY"

        # hotstart
        if self.restart is not None:
            hotout = int((self.restart - self.rdate).total_seconds() / (params["tstep"] * params["tstep_graph"]))
            params["hotstart"] = True
            params["restart_tstep"] = hotout

        if self.monitor:
            params["monitor"] = True

        # custom fortran file
        if self.fortran:
            params["fortran"] = True
            if os.path.isfile(self.fortran):
                file = os.path.split(self.fortran)[1]
                os.makedirs(os.path.join(path, "user_fortran"), exist_ok=True)
                shutil.copy(self.fortran, os.path.join(path, "user_fortran", file))
            elif os.path.isdir(self.fortran):
                files_to_copy = [os.path.basename(x) for x in glob.glob(f"{self.fortran}/*.f*")]
                files_to_copy += [os.path.basename(x) for x in glob.glob(f"{self.fortran}/*.F*")]
                dest_name = os.path.join(path, "user_fortran")
                copy_files(dest_name, self.fortran, files_to_copy)
            else:
                raise ValueError(f"Couldn't find valid FORTRAN files in {self.fortran}")
            logger.info(f"Copied FORTRAN files to {path}")

        # update
        if dic:
            for key in dic.keys():
                params[key] = dic[key]
        if "params" in kwargs.keys():
            for key in kwargs["params"].keys():
                params[key] = kwargs["params"][key]

        self.params = params

        if output:
            # save params
            # ---------------------------------------------------------------------
            logger.info("output " + self.module + " CAS file ...\n")
            # ---------------------------------------------------------------------
            write_cas(path, self.module, params)

    # ============================================================================================
    # METEO
    # ============================================================================================
    def force(self, **kwargs):
        meteo_source = get_value(self, kwargs, "meteo_source", None)

        kwargs.update({"meteo_source": meteo_source})

        flag = get_value(self, kwargs, "update", ["all"])

        z = {**self.__dict__, **kwargs}  # merge self and possible kwargs

        # check if files exist
        if flag:
            if ("meteo" in flag) | ("all" in flag):
                self.meteo = pmeteo.Meteo(**z)
            else:
                logger.info("skipping meteo ..\n")
        else:
            self.meteo = pmeteo.Meteo(**z)

    def to_force(self, geo, outpath):
        # # WRITE METEO FILE
        logger.info("saving meteo file.. ")
        # define file names
        atm_zarr = os.path.join(outpath, "atm.zarr")
        geo_mesh = xr.open_dataset(geo, engine="selafin")
        meteoFile = os.path.join(outpath, "input_wind.slf")
        # temp zarr file
        write_meteo_on_mesh(
            self.meteo.Dataset,
            geo_mesh,
            atm_zarr,
            50,
            len(geo_mesh.x),
            gtype=self.gtype,
            ttype=self.ttype,
            input360=self.input360,
        )
        # zarr to selafin
        write_meteo_selafin(meteoFile, atm_zarr)

    # ============================================================================================
    # TPXO
    # ============================================================================================
    def tpxo(self, tpxo_suffix="tpxo9_atlas_30_v4", **kwargs):
        tpxo_source = get_value(self, kwargs, "tpxo_source", None)
        path = get_value(self, kwargs, "rpath", "./telemac")

        lat_min = np.round(self.lat_min - 1 if self.lat_min > -89 else self.lat_min, 2)
        lat_max = np.round(self.lat_max + 1 if self.lat_max < 89 else self.lat_max, 2)
        lon_min = np.round(self.lon_min - 1, 2)
        lon_max = np.round(self.lon_max + 1, 2)

        if lon_min < 0 and lon_max < 0:
            lon_min += 360
            lon_max += 360

        if tpxo_source is None:
            # ---------------------------------------------------------------------
            logger.error("model not set properly: define tpxo_source\n")
            raise ValueError("model not set properly: define tpxo_source")
            # ---------------------------------------------------------------------
        os.makedirs(os.path.join(path, "TPXO"), exist_ok=True)

        logger.info("extracting tides..\n")
        #
        lines = []
        lines.append(f"{path}/TPXO/input_sources !")
        lines.append(f"{path}/TPXO/Model_Local        ! Local area Model name")
        lines.append(f"{lat_min} {lat_max}     ! Lat limits (degrees N)")
        lines.append(f"{lon_min} {lon_max}     ! Lon limits (degrees E, range -180 +360)")
        with open(f"{path}/TPXO/setup.local", "w") as file:
            file.write("\n".join(lines))

        # input sources
        lines = []
        lines.append(f"{tpxo_source}/h_*_{tpxo_suffix}")
        lines.append(f"{tpxo_source}/u_*_{tpxo_suffix}")
        lines.append(f"{tpxo_source}/grid_{tpxo_suffix}")
        with open(f"{path}/TPXO/input_sources", "w") as file:
            file.write("\n".join(lines))

        # output files
        lines = []
        lines.append(f"{path}/TPXO/h_LOCAL")
        lines.append(f"{path}/TPXO/uv_LOCAL")
        lines.append(f"{path}/TPXO/grid_LOCAL")
        with open(f"{path}/TPXO/Model_Local", "w") as file:
            file.write("\n".join(lines))

        os.system(f"cp {DATA_PATH}/extract_local_model {path}/TPXO/extract_local_model")
        os.system(f"{path}/TPXO/extract_local_model<{path}/TPXO/setup.local")

    # ============================================================================================
    # DEM
    # ============================================================================================
    def bath(self, **kwargs):
        #       z = self.__dict__.copy()

        if self.mesh.Dataset is not None:
            kwargs["grid_x"] = self.mesh.Dataset.SCHISM_hgrid_node_x.values
            kwargs["grid_y"] = self.mesh.Dataset.SCHISM_hgrid_node_y.values

        dpath = get_value(self, kwargs, "dem_source", None)

        kwargs.update({"dem_source": dpath})

        flag = get_value(self, kwargs, "update", ["all"])
        # check if files exist
        if flag:
            if ("dem" in flag) | ("all" in flag):
                kwargs.update(
                    {
                        "lon_min": self.lon_min,
                        "lat_min": self.lat_min,
                        "lon_max": self.lon_max,
                        "lat_max": self.lat_max,
                    }
                )
                self.dem.Dataset = pdem.dem_on_mesh(self.dem.Dataset, **kwargs)

                if kwargs.get("adjust_dem", True):
                    coastline = kwargs.get("coastlines", None)
                    if coastline is None:
                        logger.warning("coastlines not present, aborting adjusting dem\n")
                    elif "adjusted" in self.dem.Dataset.data_vars.keys():
                        logger.info("Dem already adjusted\n")
                    elif "fval" in self.dem.Dataset.data_vars.keys():
                        logger.info("Dem already adjusted\n")
                    else:
                        self.dem.adjust(coastline, **kwargs)
                try:
                    self.get_depth_from_dem()
                except AttributeError as e:
                    logger.info("Keeping bathymetry in hgrid.gr3 due to {}\n".format(e))
            else:
                logger.info("dem from mesh file\n")

    def get_depth_from_dem(self):
        try:
            logger.info("Dem already adjusted\n")
            bat = -self.dem.Dataset.fval.values.astype(float)  # minus for hydro run
            if np.isinf(bat).sum() != 0:
                raise Exception("Bathymetry contains Infs")
            if np.isnan(bat).sum() != 0:
                raise Exception("Bathymetry contains NaNs")

        except:
            bat = -self.dem.Dataset.ival.values.astype(float)  # minus for hydro run

        self.mesh.Dataset.depth.loc[: bat.size] = bat
        logger.info("updating bathymetry ..\n")

    @staticmethod
    def mesh_to_slf(
        x, y, z, tri, outpath, tag="telemac2d", chezy=None, manning=0.027, friction_type="chezy", **kwargs
    ):
        #
        ds = xr.Dataset(
            {
                "B": (("time", "node"), z[np.newaxis, :]),
                # Add other variables as needed
            },
            coords={
                "x": ("node", x),
                "y": ("node", y),
                "time": [pd.Timestamp.now()],
            },
        )
        ds.attrs["ikle2"] = tri + 1
        # bottom friction only in the case of TELEMAC2D
        if tag == "telemac2d":
            if chezy:
                c = np.ones(len(z)) * chezy
            else:
                if friction_type == "chezy":
                    c = (abs(z) ** (1 / 6)) / manning
                else:
                    print("only Chezy implemented so far! ")
                    sys.exit()
            ds["W"] = xr.Variable(dims=["node"], data=c)
            logger.info("Manning file created..\n")
        ds.selafin.write(outpath)

    def to_slf(self, outpath, global_=True, flip_=True, friction_type="chezy", **kwargs):
        corrections = get_value(self, kwargs, "mesh_corrections", {"reverse": [], "remove": []})

        X = self.mesh.Dataset.SCHISM_hgrid_node_x.data
        Y = self.mesh.Dataset.SCHISM_hgrid_node_y.data
        # depth
        Z = -self.mesh.Dataset.depth.data
        nan_mask = np.isnan(Z)
        inf_mask = np.isinf(Z)
        if np.any(nan_mask) or np.any(inf_mask):
            Z[nan_mask] = 0
            Z[inf_mask] = 0
            logger.info("Z contains . Cleaning needed.")
        # connectivity
        IKLE2 = self.mesh.Dataset.SCHISM_hgrid_face_nodes.data
        remove(outpath)

        if global_:
            # suppress the pole triangles (for now)
            IKLE2 = remove_pole_triangles(X, Y, IKLE2, corrections)
        if flip_:
            # adjust triangles orientation on the dateline
            IKLE2 = flip_triangles(X, Y, IKLE2, corrections)

        if IKLE2.shape != self.mesh.Dataset.SCHISM_hgrid_face_nodes.data.shape:
            reassign_ikle = True
        else:
            if not np.array_equal(IKLE2, self.mesh.Dataset.SCHISM_hgrid_face_nodes.data):
                reassign_ikle = True
            else:
                reassign_ikle = False

        if reassign_ikle:
            self.mesh.Dataset["SCHISM_hgrid_face_nodes"] = xr.Variable(
                ("nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"),
                IKLE2,
            )

        # write mesh
        chezy = get_value(self, kwargs, "chezy", None)
        manning = get_value(self, kwargs, "manning", 0.027)
        self.mesh_to_slf(X, Y, Z, IKLE2, outpath, self.module, chezy, manning, friction_type)

    # ============================================================================================
    # EXECUTION
    # ============================================================================================
    def create(self, **kwargs):
        if not kwargs:
            kwargs = self.__dict__.copy()

        # Set background dem as scale for mesh generation
        dpath = get_value(self, kwargs, "dem_source", None)
        if dpath:
            self.dem = pdem.Dem(**kwargs)
            # kwargs.update({"dem_source": self.dem.Dataset})
        else:
            logger.info("no dem available\n")

        # Mesh
        self.mesh = pmesh.set(type="tri2d", **kwargs)

        # set lat/lon from file
        if self.mesh.Dataset is not None:
            kwargs.update({"lon_min": self.mesh.Dataset.SCHISM_hgrid_node_x.values.min()})
            kwargs.update({"lon_max": self.mesh.Dataset.SCHISM_hgrid_node_x.values.max()})
            kwargs.update({"lat_min": self.mesh.Dataset.SCHISM_hgrid_node_y.values.min()})
            kwargs.update({"lat_max": self.mesh.Dataset.SCHISM_hgrid_node_y.values.max()})

            self.lon_min = self.mesh.Dataset.SCHISM_hgrid_node_x.values.min()
            self.lon_max = self.mesh.Dataset.SCHISM_hgrid_node_x.values.max()
            self.lat_min = self.mesh.Dataset.SCHISM_hgrid_node_y.values.min()
            self.lat_max = self.mesh.Dataset.SCHISM_hgrid_node_y.values.max()

        # get bathymetry
        if len(self.mesh.Dataset.depth) == len(self.mesh.Dataset.SCHISM_hgrid_node_y):
            if np.all(self.mesh.Dataset.depth == 0):
                logger.info("bathy is set to 0, interpolating dem on mesh")
                self.bath(**kwargs)
            logger.info("found bathy in mesh file\n")
        elif self.mesh.Dataset is not None:
            self.bath(**kwargs)
        else:
            raise ValueError("No bathymetry found")

        # get meteo
        if self.atm:
            self.force(**kwargs)

        # get tide
        if self.tide:
            self.tpxo(**kwargs)

        self.config(**kwargs)

    def output(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./telemac/")
        global_ = get_value(self, kwargs, "global", True)
        chezy = self.params.get("chezy", None)
        flip_ = get_value(self, kwargs, "flip", True)

        if not os.path.exists(path):
            os.makedirs(path)

        # Mesh related files
        if self.mesh.Dataset is not None:
            try:
                self.get_depth_from_dem()
            except AttributeError as e:
                logger.info("Keeping bathymetry from hgrid.gr3 ..\n")

            logger.info("saving geometry file.. ")
            geo = os.path.join(path, "geo.slf")
            cli = os.path.join(path, "geo.cli")
            self.to_slf(geo, global_=global_, flip_=flip_, chezy=chezy)
            self.mesh.to_file(filename=os.path.join(path, "hgrid.gr3"))
            write_netcdf(self.mesh.Dataset, geo)

            # WRITE METEO FILE
            logger.info("saving meteo file.. ")
            meteoFile = os.path.join(path, "input_wind.slf")
            atm_zarr = os.path.join(path, "atm.zarr")
            geo_mesh = xr.open_dataset(geo, engine="selafin")
            if isinstance(self.meteo, list):
                self.meteo = pmeteo.Meteo(self.meteo_source)
            else:
                pass

            if self.meteo:
                write_meteo_on_mesh(
                    self.meteo.Dataset,
                    geo_mesh,
                    atm_zarr,
                    50,
                    len(geo_mesh.x),
                    gtype=self.gtype,
                    ttype=self.ttype,
                    input360=self.input360,
                )
                write_meteo_selafin(meteoFile, atm_zarr)

            # WRITE BOUNDARY FILE
            logger.info("saving boundary file.. ")
            domain = export_cli(self.mesh.Dataset, geo, cli, tel_module=self.module)
            self.params["N_bc"] = len(domain[0])

            # WRITE CAS FILE
            logger.info("saving CAS file.. ")
            write_cas(path, self.module, self.params)

        # ---------------------------------------------------------------------
        logger.info("output done\n")
        # ---------------------------------------------------------------------

    def run(self, api=True, **kwargs):
        calc_dir = get_value(self, kwargs, "rpath", "./telemac/")
        cpu = get_value(self, kwargs, "cpu", os.cpu_count())
        cwd = os.getcwd()

        # ---------------------------------------------------------------------
        logger.info("executing model\n")
        # ---------------------------------------------------------------------
        comm = MPI.COMM_WORLD

        cas_file = self.module + ".cas"
        os.chdir(calc_dir)
        if not tools.is_mpirun_installed():
            logger.warning("mpirun is not installed, ending.. \n")
            return

        if self.fortran:
            user_fortran = "user_fortran"
        else:
            user_fortran = None

        if api:
            if self.module == "telemac2d":
                # Creation of the instance Telemac2d
                study = Telemac2d(cas_file, user_fortran=user_fortran, comm=comm, stdout=0, recompile=True)
            elif self.module == "telemac3d":
                study = Telemac3d(cas_file, user_fortran=user_fortran, comm=comm, stdout=0, recompile=True)
            elif self.module == "tomawac":
                study = Tomawac(cas_file, user_fortran=user_fortran, comm=comm, stdout=0, recompile=True)
            else:
                raise ValueError("this module", self.module, "is not implemented yet!")

            # Testing construction of variable list
            _ = study.variables

            study.set_case()
            if self.module == "tomawac":
                study.set("MODEL.RESULTFILE", "results_2D.slf")
            # Initalization
            study.init_state_default()
            ntimesteps = study.get("MODEL.NTIMESTEPS")
            pbar = ProgressBar(ntimesteps)
            for it in range(ntimesteps):
                study.run_one_time_step()
                pbar.update(it)

            # Ending the run
            study.finalize()
            pbar.finish()
        else:
            logger.info("no api mode, running telemac manually.. \n")
            os.system(f"{self.module}.py {cas_file} --ncsize {cpu} -s")
        #
        os.chdir(cwd)

    def save(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./telemac/")

        lista = [key for key, value in self.__dict__.items() if key not in ["meteo", "dem", "mesh"]]
        dic = {k: self.__dict__.get(k, None) for k in lista}

        mesh = self.__dict__.get("mesh", None)
        if isinstance(mesh, str):
            dic.update({"mesh": mesh})
        else:
            dic.update({"mesh": mesh.__class__.__name__})

        dem = self.__dict__.get("dem", None)
        if isinstance(dem, str):
            dic.update({"dem": dem})
        elif isinstance(dem, pdem.Dem):
            dic.update({"dem_attrs": dem.Dataset.elevation.attrs})

        meteo = self.__dict__.get("meteo", None)
        if isinstance(meteo, str):
            dic.update({"meteo": meteo})
        elif isinstance(meteo, pmeteo.Meteo):
            try:
                dic.update({"meteo": [meteo.Dataset.attrs]})
            except:
                dic.update({"meteo": [x.attrs for x in meteo.Dataset]})

        coastlines = self.__dict__.get("coastlines", None)
        if coastlines is not None:
            # Save the path to the serialized coastlines - #130
            coastlines_database = os.path.join(path, "coastlines.json")
            coastlines.to_file(coastlines_database)
            dic.update({"coastlines": coastlines_database})
        else:
            dic.update({"coastlines": None})

        dic["version"] = pyposeidon.__version__

        for attr, value in dic.items():
            if isinstance(value, datetime.datetime):
                dic[attr] = value.isoformat()
            if isinstance(value, pd.Timedelta):
                dic[attr] = value.isoformat()
            if isinstance(value, pd.DataFrame):
                dic[attr] = value.to_dict()
            if isinstance(value, gp.GeoDataFrame):
                dic[attr] = value.to_json()

        filename = os.path.join(path, f"{self.module}_model.json")
        json.dump(dic, open(filename, "w"), indent=4, default=myconverter)

    def execute(self, **kwargs):
        self.create(**kwargs)
        self.output(**kwargs)
        if self.monitor:
            self.set_obs()
        self.save(**kwargs)
        self.run(**kwargs)

    def read_folder(self, rfolder, **kwargs):
        self.rpath = rfolder
        geo = glob.glob(os.path.join(rfolder, "geo.slf"))
        cli = glob.glob(os.path.join(rfolder, "/geo.cli"))
        mfiles = glob.glob(os.path.join(rfolder, "/input_wind.slf"))

        mykeys = ["start_year", "start_month", "start_day"]
        sdate = [self.params["opt"][x] for x in mykeys]  # extract date info from param.nml
        sd = "-".join(str(item) for item in sdate)  # join in str
        sd = pd.to_datetime(sd)  # convert to datestamp

        sh = [self.params["opt"][x] for x in ["start_hour", "utc_start"]]  # get time attrs
        sd = sd + pd.to_timedelta("{}H".format(sh[0] + sh[1]))  # complete start_date

        ed = sd + pd.to_timedelta("{}D".format(self.params["core"]["rnday"]))  # compute end date based on rnday

        self.start_date = sd  # set attrs
        self.end_date = ed

        load_mesh = get_value(self, kwargs, "load_mesh", False)

        if load_mesh:
            try:
                self.mesh = pmesh.set(type="tri2d", mesh_file=hfile)
            except:
                logger.warning("loading mesh failed")
                pass
        else:
            logger.warning("no mesh loaded")

        load_meteo = get_value(self, kwargs, "load_meteo", False)

        # meteo
        if load_meteo is True:
            try:
                pm = []
                for key, val in mfiles.items():
                    msource = xr.open_mfdataset(val)
                    pm.append(msource)
                if len(pm) == 1:
                    pm = pm[0]
                self.meteo = pmeteo.Meteo(meteo_source=pm)

            except:
                pm = []
                for key, val in mfiles.items():
                    ma = []

                    for ifile in mfiles:
                        g = xr.open_dataset(ifile)
                        ts = "-".join(g.time.attrs["base_date"].astype(str)[:3])
                        time_r = pd.to_datetime(ts)
                        try:
                            times = time_r + pd.to_timedelta(g.time.values, unit="D").round("H")
                            g = g.assign_coords({"time": times})
                            ma.append(g)
                        except:
                            logger.warning("loading meteo failed")
                            break
                    if ma:
                        msource = xr.merge(ma)
                        pm.append(msource)

                if len(pm) == 1:
                    pm = pm[0]

                self.meteo = pmeteo.Meteo(meteo_source=pm)

        else:
            logger.warning("no meteo loaded")

    def hotstart(self, t=None, **kwargs):
        filename2d = get_value(self, kwargs, "filename2d", "out_2D.zarr")
        ppath = get_value(self, kwargs, "ppath", self.rpath)
        out2d = os.path.join(ppath, "outputs", filename2d)
        ds = xr.open_dataset(out2d)
        t_ = ds.time.values
        it = np.where(t_ == t)[0]
        xdat = ds.isel(time=it)
        hfile = f"hotstart_{t.strftime('%Y%m%d.%H')}.nc"
        logger.info("saving hotstart file\n")
        xdat.to_netcdf(os.path.join(ppath, f"outputs/{hfile}"))

    @staticmethod
    def set_attrs(ds):
        lat_coord_standard_name = "latitude"
        lon_coord_standard_name = "longitude"
        x_units = "degrees_east"
        y_units = "degrees_north"

        # set Attrs
        ds.longitude.attrs = {
            "long_name": "node x-coordinate",
            "standard_name": lon_coord_standard_name,
            "units": x_units,
            "mesh": "TELEMAC_Selafin",
        }

        ds.latitude.attrs = {
            "long_name": "node y-coordinate",
            "standard_name": lat_coord_standard_name,
            "units": y_units,
            "mesh": "TELEMAC_Selafin",
        }

        if "depth" in ds.variables:
            ds.depth.attrs = {
                "long_name": "Bathymetry",
                "units": "meters",
                "positive": "down",
                "mesh": "TELEMAC_Selafin",
                "location": "node",
            }

        ds.time.attrs = {
            "long_name": "Time",
            "base_date": pd.Timestamp(*ds.attrs["date_start"]).strftime("%Y-%m-%d %H:%M:%S"),
            "standard_name": "time",
        }

        # Dataset Attrs
        ds.attrs = {
            "Conventions": "CF-1.0, UGRID-1.0",
            "title": "TELEMAC Model output",
            "source": "TELEMAC model output version vV8P4",
            "references": "https://gitlab.pam-retd.fr/otm/telemac-mascaret",
            "history": "created by pyposeidon",
            "comment": "TELEMAC Model output",
            "type": "TELEMAC Model output",
        }
        return ds

    @staticmethod
    def slf_to_xarray(file: str, res_type: ResultTypes = "2D"):
        """
        Open a Selafin file and convert it to an xarray dataset.

        Args:
            file (str): Path to the Selafin file.
            type (str, optional): Type of the Selafin file. Defaults to "2D".

        Returns:
            xr.Dataset: The xarray dataset.
        """
        assert res_type in ("1D", "2D", "3D")
        ds = xr.open_dataset(file, engine="selafin")
        dic = {}
        varsn = normalize_varnames(ds.variables)
        for var, varn in zip(ds.variables, varsn):
            dic.update({var: varn})
        ds = ds.rename(dic)
        ds["triface_nodes"] = xr.Variable(("triface", "three"), np.array(ds.attrs["ikle2"]) - 1)
        if res_type == "1D":
            # bug in TELEMAC coord: CONVERT BACK FROM MERCATOR (NOT FIXED YET)
            x2, y2 = spherical2longlat(ds.longitude.values, ds.latitude.values, 0, 0)
            ds["lon"] = xr.Variable(("node"), x2)
            ds["lat"] = xr.Variable(("node"), y2)
        return ds

    def results(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./telemac/")
        filename = get_value(self, kwargs, "filename", "stations.zarr")
        filename2d = get_value(self, kwargs, "filename2d", "out_2D.zarr")
        remove_zarr = get_value(self, kwargs, "remove_zarr", True)
        chunk = get_value(self, kwargs, "chunk", None)

        logger.info("get combined 2D netcdf files \n")
        # check for new IO output
        res = os.path.join(path, "results_2D.slf")
        xc = self.slf_to_xarray(res)
        xc = self.set_attrs(xc)

        os.makedirs(os.path.join(path, "outputs"), exist_ok=True)
        out2d = os.path.join(path, "outputs", filename2d)
        remove(out2d)
        export_xarray(xc, out2d, chunk=chunk, remove_dir=remove_zarr)

        if self.monitor:
            logger.info("export observations file\n")

            res_1D = os.path.join(path, "results_1D.slf")
            xc = self.slf_to_xarray(res_1D, res_type="1D")
            xc = self.set_attrs(xc)

            out_obs = os.path.join(path, "outputs", filename)
            remove(out_obs)
            export_xarray(xc, out_obs, chunk=chunk, remove_dir=remove_zarr)
        logger.info("done with output results files 1D/2D\n")

    def set_obs(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./telemac/")
        tg_database = get_value(self, kwargs, "obs", None)
        max_dist = get_value(self, kwargs, "max_dist", np.inf)
        id_str = get_value(self, kwargs, "id_str", "ioc_code")
        offset = get_value(self, kwargs, "offset", 0)

        if tg_database:
            logger.info("get stations from {}\n".format(tg_database))
            tg = pd.read_csv(tg_database, index_col=0).reset_index(drop=True)
        else:
            logger.info("get stations using searvey\n")
            geometry = get_value(self, kwargs, "geometry", None)
            if geometry == "global":
                tg = ioc.get_ioc_stations()
            else:
                geo_box = shapely.geometry.box(self.lon_min, self.lat_min, self.lon_max, self.lat_max)
                tg = ioc.get_ioc_stations(region=geo_box)
            tg = tg.reset_index(drop=True)
            ### save in compatible to searvey format
            tg["country"] = tg.country.values.astype("str")  # fix an issue with searvey see #43 therein
            logger.info("save station DataFrame \n")
            tg_database = os.path.join(path, "stations.json")
            tg.to_file(tg_database)
        ##### normalize to be used inside pyposeidon
        tgn = normalize_column_names(tg.copy())

        ##### make sure lat/lon are floats
        tgn = tgn.astype({"latitude": float, "longitude": float, "location": str})
        ## FOR TIDE GAUGE MONITORING

        logger.info("set in-situ measurements locations \n")

        if self.mesh == "tri2d":
            self.mesh = pmesh.set(type="tri2d", mesh_file=self.mesh_file)

        if not self.mesh.Dataset:
            logger.warning("no mesh available skipping \n")
            return

        mesh = pd.DataFrame(
            np.array([self.mesh.Dataset.SCHISM_hgrid_node_x.values, self.mesh.Dataset.SCHISM_hgrid_node_y.values]).T,
            columns=["lon", "lat"],
        )
        points = pd.DataFrame(np.array([tgn.longitude.values, tgn.latitude.values]).T, columns=["lon", "lat"])
        df = find_nearest_nodes(mesh, points, 1)
        df = df[df.distance < max_dist]

        if len(df) == 0:
            logger.warning("no observations available\n")

        df["z"] = 0
        df.index += 1
        df["ind"] = df.index
        df["unique_id"] = tgn[id_str].values
        if id_str != "unique_id":
            df[id_str] = tgn[id_str].values
        df["depth"] = self.mesh.Dataset.depth[df.mesh_index.values]

        # convert to MERCATOR coordinates
        # dirty fix (this needs to be fixed in TELEMAC directly)
        x, y = longlat2spherical(df["lon"], df["lat"], 0, 0)
        df["x"] = x
        df["y"] = y

        logger.info("write out stations.csv file \n")
        df.to_csv(os.path.join(path, "stations.csv"), index=False)

        logger.info("write out stations.in file \n")
        # output to file
        serialize_stations(
            df,
            os.path.join(path, "station.in"),
            format="telemac",
            duration=self.params["duration"],
            timestep=self.params["tstep"],
            offset=offset,
        )
        self.stations = os.path.join(path, "stations.csv")

    def get_output_data(self, **kwargs):
        dic = self.__dict__

        dic.update(kwargs)

        self.data = data.get_output(**dic)

    def get_station_obs_data(self, **kwargs):
        self.station_obs_data = pd.read_csv(self.stations)

    def open_thalassa(self, **kwargs):
        # open a Thalassa instance to visualize the output
        return
