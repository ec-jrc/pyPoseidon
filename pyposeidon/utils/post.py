import os
import pandas as pd
import xarray as xr
import numpy as np

import pyposeidon
from pyposeidon.utils import data
from pyposeidon.utils.vals import obs
from pyposeidon.utils.statistics import get_stats

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
        st = b.data.time_series
        b.data.Dataset.to_netcdf(rpath + "/{}.nc".format(b.start_date.strftime("%Y%m%d%H")))  # save netcdf

        h = []
        for l in range(freq):
            from_date = b.start_date + pd.to_timedelta("{}H".format(l * 12)) + pd.to_timedelta("1S")
            to_date = b.start_date + pd.to_timedelta("{}H".format((l + 1) * 12))
            h.append(st.sel(time=slice(from_date, to_date), drop=True))

        dr = xr.concat(h, dim="lead")

        c.append(dr)

    cc = xr.concat(c, dim="time")

    # Observation/Stats
    logger.info("Retrieve observation data for station points\n")
    odata = obs(**b.__dict__, sa_date=b.date, se_date=b.end_date)

    logger.info("Compute general statistics for station points\n")
    tgs = []
    sts = []
    for idx, val in odata.locations.iterrows():

        obs_ = odata.iloc(idx)  # Get observational data

        sim = st.isel(node=idx).elev.to_dataframe().drop("node", axis=1)

        stable = get_stats(sim, obs_)  # Do general statitics

        tgs.append(obs_.Surge)  # store
        sts.append(stable)

    logger.info("Construct station data output\n")
    # Construct Thalassa file
    stations = pd.DataFrame(b.stations)
    stations = stations.reset_index(drop=True)
    stations.index.name = "node"
    stations["id"] = odata.locations.ID.values
    stations["group"] = odata.locations.Group.values
    stations["longitude"] = odata.locations.longitude.values
    stations["latitude"] = odata.locations.latitude.values

    sd = b.data.time_series.rename({"elev": "elevation"})

    oo = pd.concat(tgs, keys=np.arange(len(tgs)), names=["node", "time"])
    oo.name = "observation"

    stats = pd.DataFrame(sts)
    stats.index.name = "node"

    cf = cc.rename({"elev": "forecast", "time": "ftime"})

    vdata = xr.merge([oo.to_xarray(), stats.to_xarray(), stations.to_xarray(), sd, cf])

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
