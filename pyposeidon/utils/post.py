import os
import sys
import pandas as pd
import geopandas as gp
import xarray as xr
import numpy as np
from tqdm.auto import tqdm

import pyposeidon
from pyposeidon.utils import data
from pyposeidon.utils.statistics import get_stats
from pyposeidon.utils.obs import get_obs_data
from pyposeidon.utils.seam import to_2d

# from pyposeidon.utils.detide import get_ss

import logging

logger = logging.getLogger(__name__)


def get_encoding(ename):
    return {
        ename: {"zlib": True, "complevel": 1},
    }


def get_bogus(temp, idfs, ntime):
    return xr.DataArray(
        temp,
        coords={
            "node": (["node"], [idfs]),
            "time": ntime,
        },
    )


def compute_obs(stations, start_time, end_time):
    logger.info("Retrieve observation data for station points\n")
    odata = get_obs_data(stations=stations, start_time=start_time, end_time=end_time)

    dfs = [x for x in stations.ioc_code if x not in odata.ioc_code]
    idx = stations[stations["ioc_code"].isin(dfs)].index

    logger.info("Normalize observation data for station points\n")

    od = []
    for inode in tqdm(odata.ioc_code.values):
        oi = odata.sel(ioc_code=inode)

        for var in oi.data_vars:
            if oi[var].isnull().all().values == True:
                oi = oi.drop(var)

        var = [k for k, v in oi.data_vars.items() if v.dims == ("time",)][0]

        obs = oi[var].to_dataframe().drop(["ioc_code"], axis=1)  # Get observational data

        # de-tide obs
        #        if not obs[var].dropna().empty | (obs[var].dropna() == obs[var].dropna()[0]).all():
        #            obs = get_ss(obs, oi.lat.values)
        #            oi[var].values = obs.elev.values

        od.append(oi[var].rename("elev_obs"))

    ods = xr.concat(od, dim="ioc_code").rename({"ioc_code": "node"})  # obs

    ## add bogus observations in order to keep array shape

    ntime = odata.isel(ioc_code=0).time.data
    temp = np.array([np.NaN] * ntime.shape[0])[np.newaxis, :]
    xdfs = []
    for idfs in dfs:
        xdfs.append(get_bogus(temp, idfs, ntime))

    ods = xr.concat([ods] + xdfs, dim="node")

    return ods


def save_ods(ods, **kwargs):
    rpath = kwargs.get("rpath", "./thalassa/")

    # obs data
    logger.info("saving observations data\n")
    obs_file = os.path.join(rpath, f"searvey")
    if not os.path.exists(obs_file):
        ods.to_zarr(store=obs_file, mode="a")
    else:
        ods.to_zarr(store=obs_file, mode="a", append_dim="time")

    return


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
    b = pyposeidon.model.read(folder + "/{}_model.json".format(tag))
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

    logger.info("Saving combined netcdf file for folder {}\n".format(folder))
    output_file = os.path.join(
        rpath,
        f"{b.start_date.strftime('%Y%m%d%H')}.nc",
    )
    out.to_netcdf(output_file, encoding=encoding)  # save netcdf

    stations = gp.GeoDataFrame.from_file(b.obs)

    ## save results in lead chunks
    lpath = os.path.join(rpath, "skill")

    for l in range(freq):
        from_date = b.start_date + pd.to_timedelta("{}H".format(l * 12)) + pd.to_timedelta("1S")
        to_date = b.start_date + pd.to_timedelta("{}H".format((l + 1) * 12))
        h = st.sel(time=slice(from_date, to_date), drop=True)
        h = h.rename({"elev": "elev_fct", "time": "ftime"})
        h = h.assign_coords({"node": stations.ioc_code.to_list()})
        leadfile = os.path.join(lpath, f"lead{l}")
        if not os.path.exists(leadfile):
            h.to_zarr(store=leadfile, mode="a")
        else:
            h.to_zarr(store=leadfile, mode="a", append_dim="ftime")

    locations = stations.ioc_code.values

    # save to files
    logger.info("Construct station simulation data output\n")
    # Construct Thalassa file

    stp = stations.to_xarray().rename({"index": "node"})  # stations
    stp = stp.assign_coords({"node": locations})

    st_ = st.assign_coords(node=locations)
    st_ = st_.rename({"elev": "elev_sim", "time": "stime"})

    vdata = xr.merge([stp, st_])

    vdata = vdata.drop_vars("geometry")

    logger.info("Save station simulation data output\n")

    output_path = os.path.join(rpath, filename)

    vdata.to_netcdf(output_path)
    logger.info(f"..done with {filename} file\n")

    # get observations
    ods = compute_obs(stations, b.start_date, b.end_date)

    save_ods(ods, **kwargs)

    # compute stats
    st = st.assign({"ioc_code": ("node", stations["ioc_code"])})
    sts = compute_stats(st, ods)

    save_stats(sts, stations, **kwargs)

    logger.info(f"post processing complete for folder {folder}\n")


def compute_stats(st, ods):
    logger.info("Compute general statistics for station points\n")

    sts = []
    for inode in tqdm(st.ioc_code.values):
        isim = st.where(st.ioc_code == inode).dropna(dim="node")
        sim = isim.elev.to_dataframe().droplevel(1)

        obs = ods.sel(node=inode).to_dataframe().drop("node", axis=1)

        stable = get_stats(sim, obs)  # Do general statitics

        sts.append(stable)

    return sts


def save_stats(sts, stations, **kwargs):
    rpath = kwargs.get("rpath", "./thalassa/")

    logger.info("Save stats\n")

    locations = stations.ioc_code.values

    stats = pd.DataFrame(sts)
    stats.index.name = "node"
    stats_ = stats.to_xarray().assign_coords({"node": locations})  # stats

    stp = stations.to_xarray().rename({"index": "node"})  # stations
    stp = stp.assign_coords({"node": locations})

    sdata = xr.merge([stp, stats_])

    sdata = sdata.drop_vars("geometry")

    rpath = kwargs.get("rpath", "./thalassa/")

    output_path = os.path.join(rpath, "stats.nc")

    sdata.to_netcdf(output_path)
    logger.info(f"..done with stats file\n")
