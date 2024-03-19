import os
import sys
import pandas as pd
import geopandas as gp
import xarray as xr
import numpy as np
from tqdm.auto import tqdm
import glob
import numcodecs
import tarfile
import shutil

import pyposeidon
from pyposeidon.tools import to_geodataframe
from pyposeidon.utils import data
from pyposeidon.utils.statistics import get_stats
from pyposeidon.utils.obs import get_obs_data
from pyposeidon.utils.seam import to_2d

# from pyposeidon.utils.detide import get_ss

valid_sensors = [
    "rad",
    "prs",
    "enc",
    "pr1",
    "PR2",
    "pr2",
    "ra2",
    "bwl",
    "wls",
    "aqu",
    "ras",
    "pwl",
    "bub",
    "enb",
    "atm",
    "flt",
    "ecs",
    "stp",
    "prte",
    "prt",
    "ra3",
]

import logging

logger = logging.getLogger(__name__)


def get_encoding(ename):
    return {
        ename: {"zlib": True, "complevel": 1},
    }


def remove(path):
    try:
        if os.path.isfile(path):
            os.remove(path)  # Remove a file
        elif os.path.isdir(path):
            if not os.listdir(path):  # Check if the directory is empty
                os.rmdir(path)
            else:
                shutil.rmtree(path)
    except OSError as e:
        print(f"Error: {e.strerror}")


def tar_directory(input_directory, output_filename, compress=False):
    """
    Tar a Zarr directory.

    Parameters
    ----------
    input_directory : str
        The path to the Zarr directory to archive.
    output_filename : str
        The path and filename for the output archive.
    compress : bool, optional
        Whether to compress the archive with gzip. Default is True.
    """
    mode = "w:gz" if compress else "w"
    with tarfile.open(output_filename, mode) as tar:
        tar.add(input_directory, arcname=os.path.basename(input_directory))


def export_xarray(ds, filename_out, chunk=None, remove_dir=False):
    """
    Export an xarray dataset to netcdf or zarr format.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray dataset to export.
    filename_out : str
        The path and filename of the output file.
    chunk : dict, optional
        The chunk size to use when saving the dataset to zarr format, by default None.
        example: chunk = {'nodes': 1000, 'time': 100}

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the filename does not have a .nc or .zarr extension.
    """
    if not filename_out.endswith(".nc") and not filename_out.endswith(".zarr"):
        raise ValueError("Invalid filename extension. Must be .nc or .zarr")
    if filename_out.endswith(".nc"):
        encoding = {}
        default_encoding = {"zlib": True, "complevel": 1}
        for varname in list(ds.variables):
            encoding[varname] = default_encoding
            if chunk:
                chunksizes = dict(ds[varname].sizes)
                for key in chunk:
                    if key in chunksizes:
                        chunksizes[key] = chunk[key]
                encoding[varname]["chunksizes"] = list(chunksizes.values())
        ds.to_netcdf(filename_out, encoding=encoding)
    elif filename_out.endswith(".zarr"):
        compressor = numcodecs.Blosc(cname="zstd", clevel=1)
        if chunk is not None:
            ds = ds.chunk(chunk)
        encoding = {}
        for varname in list(ds.variables):
            encoding.update({varname: {"compressor": compressor}})
        ds.to_zarr(filename_out, encoding=encoding, consolidated=True)
        tar_directory(filename_out, filename_out + ".tar")
        if remove_dir:
            remove(filename_out)
            logger.info(f"Removed directory {filename_out}")


def save_leads(stations, st, start_date, dt, leads, rpath="./skill/"):
    ## save results in lead chunks

    for l in range(leads):
        from_date = start_date + pd.to_timedelta("{}H".format(l * dt)) + pd.to_timedelta("1S")
        to_date = start_date + pd.to_timedelta("{}H".format((l + 1) * dt))
        h = st.sel(time=slice(from_date, to_date), drop=True)
        h = h.rename({"elev": "elev_fct", "time": "ftime"})
        leadfile = os.path.join(rpath, f"lead{l}")
        if not os.path.exists(leadfile):
            h.to_zarr(store=leadfile, mode="a")
        else:
            h.to_zarr(store=leadfile, mode="a", append_dim="ftime")

    return


def gather_obs_data(stations, start_time, end_time, rpath="./thalassa/obs/"):
    logger.info(f"Retrieve observation data for station points from {start_time} to {end_time}\n")
    odata = get_obs_data(stations=stations, start_time=start_time, end_time=end_time)

    logger.info("Save observation data for station locations\n")

    for inode in tqdm(odata.id.values):
        oi = odata.sel(id=inode)

        var = [k for k, v in oi.data_vars.items() if "time" in v.dims]

        df = oi[var].to_dataframe().drop("id", axis=1)
        df_ = df.dropna(axis=1, how="all")  # drop all nan columns

        file_path = os.path.join(rpath, f"{inode}.parquet")

        if os.path.isfile(file_path):
            obs = pd.read_parquet(file_path, engine="fastparquet")
            # make sure there is output
            if df_.empty:
                cols = [x for x in var if x in obs.columns]
                df = df[cols]
                df.to_parquet(file_path, engine="fastparquet", append=True)
            else:
                df = df_
                out = pd.concat([obs, df]).dropna(how="all")  # merge
                out.to_parquet(file_path, engine="fastparquet")
        else:
            # make sure there is output
            if df_.empty:
                df = df[[df.columns[0]]]
            else:
                df = df_

            df.to_parquet(file_path, engine="fastparquet")

    return


def to_stats(st, rpath="./thalassa/", opath="./thalassa/obs/"):
    logger.info("compute general statistics for station points\n")

    ids = st.id.values

    sts = []
    for inode in tqdm(ids):
        sim = st.sel(id=inode).elev_sim.to_dataframe().drop("id", axis=1)

        filename = os.path.join(opath, f"{inode}.parquet")
        obs = pd.read_parquet(filename, engine="fastparquet")
        obs_ = obs.dropna(axis=1, how="all")  # drop all nan columns

        if obs_.empty:
            obs = obs[[obs.columns[0]]]
        else:
            cols = [x for x in obs_.columns if x in valid_sensors]
            obs = obs_[[cols[0]]]  # just choose one for now

        if obs.dropna().empty:
            logger.warning(f"Observation data not available for {inode} station")

        stable = get_stats(sim, obs)  # Do general statitics

        sts.append(stable)

    logger.info("save stats\n")

    stats = pd.DataFrame(sts)
    stats.index.name = "id"
    stats = stats.to_xarray().assign_coords({"id": ids})  # stats

    output_path = os.path.join(rpath, "stats.nc")

    stats.to_netcdf(output_path)
    logger.info(f"..done with stats file\n")


def to_thalassa(folder, **kwargs):
    # Retrieve data
    tag = kwargs.get("tag", "schism")
    leads = kwargs.get("leads", None)
    rpath = kwargs.get("rpath", "./thalassa/")
    gglobal = kwargs.get("gglobal", False)
    to2d = kwargs.get("to2d", None)
    filename = kwargs.get("filename", "stations.nc")

    # generalize
    x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
    y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
    tes_var = kwargs.get("e", "SCHISM_hgrid_face_nodes")

    # sanity check
    if gglobal and not to2d:
        raise ValueError("global dataset but no to2d matrix is given")

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    # Get simulation data
    json_file = os.path.join(folder, "{}_model.json".format(tag))
    b = pyposeidon.model.read(json_file)
    b.get_output_data()
    b.get_station_sim_data()
    st = b.station_sim_data
    out = b.data.Dataset
    ename = [x for x in out.data_vars if x in ["elev", "elevation"]][0]

    rvars = kwargs.get("vars", [ename, "depth"])

    if gglobal:
        logger.info("converting to 2D\n")
        # Convert to 2D
        [xn, yn, tri3n] = np.load(to2d, allow_pickle=True)
        out = to_2d(out, data_vars=rvars, mesh=[xn, yn, tri3n])  # elevation

    else:
        rvars_ = rvars + [x_var, y_var, tes_var]
        out = out[rvars_]

    # Add max elevation variable
    out = out.assign(max_elev=out[ename].max("time"))
    rvars = rvars + ["max_elev"]

    # set enconding
    encoding = {}
    for var in rvars:
        encoding.update(get_encoding(var))

    logger.info("saving combined netcdf file for folder {}\n".format(folder))
    output_file = os.path.join(
        rpath,
        f"{b.start_date.strftime('%Y%m%d%H')}.nc",
    )
    out.to_netcdf(output_file, encoding=encoding)  # save netcdf

    stations = to_geodataframe(b.obs)

    # assign unique id
    if "id" not in stations.columns:
        stations["id"] = [f"IOC-{x}" for x in stations.ioc_code]
    st = st.assign_coords({"id": ("node", stations.id.values)}).swap_dims({"node": "id"}).reset_coords("node")

    logger.info("save time series depending on lead time\n")

    if leads:
        skill_path = os.path.join(rpath, "skill")
        total_hours = pd.to_timedelta(b.time_frame) / pd.Timedelta(hours=1)
        dt = total_hours / leads

        if dt % 1 == 0.0:
            save_leads(stations, st, b.start_date, dt, leads, rpath=skill_path)

        else:
            logger.warning("leads value not correct, aborting\n")

    # save sim data

    ids = stations.id.values

    stp = stations.to_xarray().rename({"index": "node"})  # stations
    stp = stp.assign_coords({"id": ("node", ids)}).swap_dims({"node": "id"}).reset_coords("node")

    st_ = st.rename({"elev": "elev_sim", "time": "stime"})

    vdata = xr.merge([stp, st_])

    vdata = vdata.drop_vars("geometry")

    logger.info("save station simulation data output\n")

    output_path = os.path.join(rpath, filename)

    # output sim data
    vdata[["node", "ioc_code", "lat", "lon", "location", "elev_sim"]].to_netcdf(output_path)
    logger.info(f"..done with {filename} file\n")


def to_obs(folder, **kwargs):
    tag = kwargs.get("tag", "schism")
    rpath = kwargs.get("rpath", "./thalassa/")

    json_file = os.path.join(folder, "{}_model.json".format(tag))
    b = pyposeidon.model.read(json_file)

    stations = to_geodataframe(b.obs)

    # assign unique id
    if "id" not in stations.columns:
        stations["id"] = [f"IOC-{x}" for x in stations.ioc_code]

    # get observations last timestamp
    obs_files_path = os.path.join(rpath, "obs/")
    if not os.path.exists(obs_files_path):
        os.makedirs(obs_files_path)

    obs_files = glob.glob(obs_files_path + "*.parquet")

    if not obs_files:
        start_date = b.start_date
    else:
        obs_ = pd.read_parquet(obs_files[0], engine="fastparquet")
        start_date = obs_.index[-1] + pd.Timedelta("60S")

    gather_obs_data(stations, start_date, b.end_date, rpath=obs_files_path)

    return
