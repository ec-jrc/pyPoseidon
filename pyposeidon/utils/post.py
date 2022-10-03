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
import pyposeidon.utils.pplot
from pyposeidon.utils.detide import get_ss

import logging

logger = logging.getLogger(__name__)


encoding = dict(
    elev={"dtype": "float32", "zlib": True, "complevel": 1, "_FillValue": -999},
    time={"dtype": "float32", "zlib": True, "complevel": 1, "_FillValue": -999},
    depth={"dtype": "float32", "zlib": True, "complevel": 1, "_FillValue": -999},
    SCHISM_hgrid_node_x={"dtype": "float32", "zlib": True, "complevel": 1, "_FillValue": -999},
    SCHISM_hgrid_node_y={"dtype": "float32", "zlib": True, "complevel": 1, "_FillValue": -999},
    SCHISM_hgrid_face_nodes={"dtype": "int32", "zlib": True, "complevel": 1, "_FillValue": 32767},
)


def to_thalassa(folders, freq=None, **kwargs):

    # Retrieve data
    tag = kwargs.get("tag", "schism")
    rpath = kwargs.get("rpath", "./thalassa/")
    gglobal = kwargs.get("gglobal", False)
    to2d = kwargs.get("to2d", None)

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    # Get simulation data
    c = []

    for folder in folders:

        b = pyposeidon.model.read(folder + "/{}_model.json".format(tag))
        b.get_output_data()
        b.get_station_sim_data()
        st = b.station_sim_data
        out = b.data.Dataset

        if gglobal:
            # Convert to 2D
            if to2d:
                [xn, yn, tri3n] = np.load(to2d, allow_pickle=True)
            else:
                logger.error("to2d conversion not given, aborting ..\n")
                sys.exit()
            sl = to_2d(out, var="elev", mesh=[xn, yn, tri3n])  # elevation
            sd = to_2d(out, var="depth", mesh=[xn, yn, tri3n])  # depth
            out = xr.merge([sl, sd])

        else:

            out = out[["elev", "depth", "SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "SCHISM_hgrid_face_nodes"]]

        # save to mp4
        logger.info("Saving animation of SSH for folder {}\n".format(folder))
        se = out.pplot.frames(var="elev", title="SSH")
        se.save(rpath + "/{}.mp4".format(b.start_date.strftime("%Y%m%d%H")))  # save mp4'

        logger.info("Saving combined netcdf file for folder {}\n".format(folder))
        out.to_netcdf(rpath + "/{}.nc".format(b.start_date.strftime("%Y%m%d%H")), encoding=encoding)  # save netcdf

        h = []
        for l in range(freq):
            from_date = b.start_date + pd.to_timedelta("{}H".format(l * 12)) + pd.to_timedelta("1S")
            to_date = b.start_date + pd.to_timedelta("{}H".format((l + 1) * 12))
            h.append(st.sel(time=slice(from_date, to_date), drop=True))

        dr = xr.concat(h, dim="lead")

        c.append(dr)

    cc = xr.concat(c, dim="time")

    cf = cc.rename({"elev": "elev_fct", "time": "ftime"})

    # Observation/Stats
    logger.info("Retrieve observation data for station points\n")
    stations = gp.GeoDataFrame.from_file(b.obs)
    odata = get_obs_data(stations=stations, start_time=b.rdate, end_time=b.end_date)

    st = st.assign({"ioc_code": ("node", stations["ioc_code"])})
    cf = cf.assign({"ioc_code": ("node", stations["ioc_code"])})

    dfs = [x for x in stations.ioc_code if x not in odata.ioc_code]
    idx = stations[stations["ioc_code"].isin(dfs)].index

    logger.info("Compute general statistics for station points\n")
    sts = []
    od = []
    for inode in tqdm(odata.ioc_code.values):

        oi = odata.sel(ioc_code=inode)

        for var in oi.data_vars:
            if oi[var].isnull().all().values == True:
                oi = oi.drop(var)

        var = [k for k, v in oi.data_vars.items() if v.dims == ("time",)][0]

        obs = oi[var].to_dataframe().drop(["ioc_code"], axis=1)  # Get observational data

        # de-tide obs
        if not obs[var].dropna().empty | (obs[var].dropna() == obs[var].dropna()[0]).all():
            obs = get_ss(obs, oi.lat.values)
            oi[var].values = obs.elev.values

        isim = st.where(st.ioc_code == inode).dropna(dim="node")
        sim = isim.elev.to_dataframe().droplevel(1)

        stable = get_stats(sim, obs)  # Do general statitics

        sts.append(stable)
        od.append(oi[var].rename("elev_obs"))

    logger.info("Construct station data output\n")
    # Construct Thalassa file

    ods = xr.concat(od, dim="ioc_code").rename({"ioc_code": "node"})  # obs

    stats = pd.DataFrame(sts)
    stats.index.name = "node"
    stats_ = stats.to_xarray().assign_coords({"node": ods.node.values})  # stats

    stations.drop(idx, inplace=True)

    stp = stations.to_xarray().rename({"index": "node"})  # stations
    stp = stp.assign_coords({"node": ods.node.values})

    st_ = st.where(st.ioc_code.isin(stats_.node), drop=True)  # sims
    st_ = st_.assign_coords(node=ods.node.values)
    st_ = st_.rename({"elev": "elev_sim", "time": "stime"})

    cf_ = cf.where(cf.ioc_code.isin(stats_.node), drop=True)  # leads
    cf_ = cf_.assign_coords(node=ods.node.values)

    sd = st.rename({"elev": "elev_sim", "time": "stime"})

    vdata = xr.merge([ods.to_dataset(), stp, stats_, cf_, st_])

    vdata = vdata.drop_vars("geometry")

    logger.info("Save station data output\n")
    vdata.to_netcdf(rpath + "/fskill.nc")
    logger.info("..done\n")


def fskill(dset, var, node=None):

    if not node:
        logger.error("Specify node\n")
        return

    lstat = []

    for l in dset.lead.values:

        obs_ = dset.sel(node=node).elev_obs  # Get observational data
        obs_ = obs_.dropna(dim="time").to_dataframe().drop("node", axis=1)

        sim = dset.sel(node=node).elev_fct.sel(lead=l).to_dataframe().drop("node", axis=1).dropna()

        stable = get_stats(sim, obs_)  # Do general statitics

        lstat.append(stable[var])

    # match values to scale
    return lstat
