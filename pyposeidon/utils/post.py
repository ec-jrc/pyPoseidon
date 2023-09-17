import os
import sys
import pandas as pd
import geopandas as gp
import xarray as xr
import numpy as np
from tqdm.auto import tqdm
import glob

import pyposeidon
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


def save_leads(stations, st, start_date, dt, freq, rpath="./skill/"):
    ## save results in lead chunks

    for l in range(freq):
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


def compute_stats(st, rpath="./thalassa/obs/"):
    logger.info("compute general statistics for station points\n")

    sts = []
    for inode in tqdm(st.id.values):
        isim = st.where(st.id == inode).dropna(dim="id")
        sim = isim.elev.to_dataframe().droplevel(1)

        filename = os.path.join(rpath, f"{inode}.parquet")
        obs = pd.read_parquet(filename, engine="fastparquet")
        obs_ = obs.dropna(axis=1, how="all")  # drop all nan columns

        if obs_.empty:
            obs = obs[[obs.columns[0]]]
        else:
            cols = [x for x in obs_.columns if x in valid_sensors]
            obs = obs_[[cols[0]]]  # just choose one for now

        stable = get_stats(sim, obs)  # Do general statitics

        sts.append(stable)

    return sts


def save_stats(sts, stations, **kwargs):
    rpath = kwargs.get("rpath", "./thalassa/")

    logger.info("save stats\n")

    ids = stations.id.values

    stats = pd.DataFrame(sts)
    stats.index.name = "id"
    stats = stats.to_xarray().assign_coords({"id": ids})  # stats

    rpath = kwargs.get("rpath", "./thalassa/")

    output_path = os.path.join(rpath, "stats.nc")

    stats.to_netcdf(output_path)
    logger.info(f"..done with stats file\n")


def to_thalassa(folder, freq=None, **kwargs):
    # Retrieve data
    tag = kwargs.get("tag", "schism")
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
        sv = []
        for var in rvars:
            isv = to_2d(out, var=var, mesh=[xn, yn, tri3n])  # elevation
            sv.append(isv)
        out = xr.merge(sv)

    else:
        rvars_ = rvars + [x_var, y_var, tes_var]
        out = out[rvars_]

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

    stations = gp.GeoDataFrame.from_file(b.obs)

    # assign unique id
    if "id" not in stations.columns:
        stations["id"] = [f"IOC-{x}" for x in stations.ioc_code]
    st = st.assign_coords({"id": ("node", stations.id.values)}).swap_dims({"node": "id"}).reset_coords("node")

    logger.info("save time series depending on lead time\n")

    skill_path = os.path.join(rpath, "skill")
    total_hours = pd.to_timedelta(b.time_frame) / pd.Timedelta(hours=1)
    dt = total_hours / freq

    if dt % 1 == 0.0:
        save_leads(stations, st, b.start_date, dt, freq, rpath=skill_path)

    else:
        logger.warning("freq not correct, aborting\n")

    # save sim data

    ids = stations.id.values

    stp = stations.to_xarray().rename({"index": "node"})  # stations
    stp = stp.assign_coords({"id": ("node", ids)}).swap_dims({"node": "id"}).reset_coords("node")

    st_ = st.rename({"elev": "elev_sim", "time": "stime"})

    vdata = xr.merge([stp, st_])

    vdata = vdata.drop_vars("geometry")

    logger.info("save station simulation data output\n")

    output_path = os.path.join(rpath, filename)

    vdata.to_netcdf(output_path)
    logger.info(f"..done with {filename} file\n")

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

    # compute stats
    logger.info("compute statistics")
    sts = compute_stats(st, rpath=obs_files_path)

    save_stats(sts, stations, **kwargs)

    logger.info(f"post processing complete for folder {folder}\n")
