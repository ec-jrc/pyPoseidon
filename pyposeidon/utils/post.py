import os
import pandas as pd
import geopandas as gp
import xarray as xr
import numpy as np

import pyposeidon
from pyposeidon.utils import data
from pyposeidon.utils.statistics import get_stats
from pyposeidon.utils.obs import get_obs_data

import logging

logger = logging.getLogger(__name__)


def to_thalassa(folders, freq=None, **kwargs):

    # Retrieve data
    tag = kwargs.get("tag", "schism")
    rpath = kwargs.get("rpath", "./thalassa/")

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    # Get simulation data
    c = []

    for folder in folders:

        b = pyposeidon.model.read(folder + "/{}_model.json".format(tag))
        b.get_output_data()
        b.get_station_sim_data()
        st = b.station_sim_data
        b.data.Dataset.to_netcdf(rpath + "/{}.nc".format(b.start_date.strftime("%Y%m%d%H")))  # save netcdf

        h = []
        for l in range(freq):
            from_date = b.start_date + pd.to_timedelta("{}H".format(l * 12)) + pd.to_timedelta("1S")
            to_date = b.start_date + pd.to_timedelta("{}H".format((l + 1) * 12))
            h.append(st.sel(time=slice(from_date, to_date), drop=True))

        dr = xr.concat(h, dim="lead")

        c.append(dr)

    cc = xr.concat(c, dim="time")

    cf = cc.rename({"elev": "forecast", "time": "ftime"})

    # Observation/Stats
    logger.info("Retrieve observation data for station points\n")
    stations = gp.GeoDataFrame(b.obs)
    odata = get_obs_data(stations=stations, start_date=b.date, end_date=b.end_date)

    logger.info("Compute general statistics for station points\n")
    sts = []
    for inode in st.node:

        obs = (
            odata.isel(ioc_code=inode).prs.to_dataframe().drop(["ioc_code", "node"], axis=1)
        )  # Get observational data

        sim = st.isel(node=inode).elev.to_dataframe().drop("node", axis=1)

        stable = get_stats(sim, obs)  # Do general statitics

        sts.append(stable)

    logger.info("Construct station data output\n")
    # Construct Thalassa file
    stats = pd.DataFrame(sts)
    stats.index.name = "node"

    stp = stations.to_xarray().rename({"index": "node"})

    od = odata.rename({"prs": "observations"}).swap_dims({"ioc_code": "node"})

    sd = st.rename({"elev": "elevation"})

    vdata = xr.merge([od, stats.to_xarray(), cf, sd])

    vdata = vdata.reset_coords("ioc_code")

    vdata = xr.merge([vdata, stp])

    vdata = vdata.drop_vars("geometry")

    logger.info("Save station data output\n")
    vdata.to_netcdf(rpath + "/fskill.nc")
    logger.info("..done\n")


def fskill(dset, var, node):

    lstat = []

    for l in dset.lead.values:

        obs_ = dset.sel(node=node).observation  # Get observational data
        obs_ = obs_.to_dataframe().drop("node", axis=1)

        sim = dset.isel(node=node).forecast.sel(lead=l).to_dataframe().drop("node", axis=1)

        stable = get_stats(sim, obs_)  # Do general statitics

        lstat.append(stable[var])

    # match values to scale
    return lstat
