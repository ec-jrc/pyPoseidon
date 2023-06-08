import datetime
import numpy as np
import pandas as pd
import utide

import logging

logger = logging.getLogger(__name__)


def get_ss(obs, lat):
    obs["seconds"] = obs.index - obs.index[0]
    obs["seconds"] = obs.seconds.dt.total_seconds().astype(int)
    obs.columns = ["elev", "seconds"]
    obs["flag"] = 0
    obs.loc[obs.elev.isna(), ["flag"]] = 2
    # interpolate and adjust bias
    bad = obs["flag"] == 2
    corrected = obs["flag"] == 1

    obs.loc[bad, "elev"] = np.nan
    obs["anomaly"] = obs["elev"] - obs["elev"].mean()
    obs["anomaly"] = obs["anomaly"].interpolate()
    #    logger.info(f"{bad.sum()} points were flagged 'bad' and interpolated")
    #    logger.info(f"{corrected.sum()} points were flagged 'corrected' and left unchanged")

    coef = utide.solve(
        obs.index,
        obs["anomaly"],
        lat=lat,
        method="ols",
        conf_int="MC",
        verbose=False,
        #    constit=["M2","K1", 'M4', 'M3', 'M8',]
    )

    tide = utide.reconstruct(obs.index, coef, verbose=False)

    dtd = obs[["elev"]].copy()
    dtd["elev"] = obs.anomaly - tide.h

    return dtd
