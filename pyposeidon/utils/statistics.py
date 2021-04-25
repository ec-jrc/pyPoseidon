#
import sys
import numpy as np
import pandas as pd
import warnings


def vtable(obsrv, model):

    if len(obsrv) != len(model):
        print("dimensions mismatch")
        return

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        #    Mean Absolute Error
        try:
            mse = np.nanmean(abs(obsrv - model))
        except:
            mse = np.NaN

        #    RMSE
        try:
            rmse = np.sqrt(np.nanmean((obsrv - model) ** 2))
        except:
            rmse = np.NaN

        #    Scatter Index
        try:
            scInd = np.sqrt(np.nanmean((obsrv - model) ** 2)) / np.nanmean(obsrv)
        except:
            scInd = np.NaN

        #   percentage RMSE
        try:
            percrmse = 100 * (rmse) / (np.nanmax(obsrv))
        except:
            percrmse = np.NaN

        #   BIAS or mean error
        # positive values indicate positively biased modelled data
        # negative values indicate negatively biased modelled data

        try:
            bias = np.nanmean(model - obsrv)
        except:
            bias = np.NaN

        # Standard deviation of residuals
        # reference: https://cirpwiki.info/wiki/Statistics

        try:
            sdr = np.sqrt(np.nanmean((model - obsrv - np.nanmean(model) + np.nanmean(obsrv)) ** 2))
        except:
            sdr = np.NaN

        # Correlation Coefficient
        # reference: https://cirpwiki.info/wiki/Statistics
        s1N = obsrv[~np.isnan(obsrv)]
        s2N = model[~np.isnan(model)]

        try:
            corcoef = np.corrcoef(s1N, s2N)[0, 1]
        except:
            corcoef = np.NaN

        # R**2
        # reference: https://cirpwiki.info/wiki/Statistics

        try:
            corcoefSq = np.corrcoef(s1N, s2N)[0, 1] ** 2
        except:
            corcoefSq = np.NaN

        # Nash-Sutcliffe Coefficient (Nash and Sutcliffe 1970)
        # Range       Qualification
        # 0.8<PS<1.0  Excellent
        # 0.6<PS<0.8  Good
        # 0.3<PS<0.6  Reasonable
        # 0<PS<0.3    Poor
        # PS<0        Bad
        # reference: https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient

        try:
            nsccoef = 1 - np.nansum((model - obsrv) ** 2) / np.nansum((obsrv - np.nanmean(obsrv)) ** 2)
        except:
            nsccoef = np.NaN

        # lamda index of Duveiller et al (2016) in ScRe
        # https://www.nature.com/articles/srep19401
        # example of application in R
        # http://eoscience.esa.int/landtraining2017/files/materials/D5P3a_I.pdf

        try:
            Xmean = np.nanmean(obsrv)
            Ymean = np.nanmean(model)
            nObs = len(obsrv)
            if corcoef >= 0:
                kappa = 0
                # print('corcoef is positive',kappa)
            else:
                kappa = 2 * abs(np.nansum((obsrv - Xmean) * (model - Ymean)))
                # print('corcoef is negative',kappa)

            Nomin = np.nansum((obsrv - model) ** 2)
            Denom = (
                np.nansum((obsrv - Xmean) ** 2)
                + np.nansum((model - Ymean) ** 2)
                + nObs * ((Xmean - Ymean) ** 2)
                + kappa
            )
            lamda_index = 1 - Nomin / Denom
        except:
            lamda_index = np.NaN

    dic = {
        "Mean Absolute Error": mse,
        "RMSE": rmse,
        "Scatter Index": scInd,
        "percentage RMSE": percrmse,
        "BIAS or mean error": bias,
        "Standard deviation of residuals": sdr,
        "Correlation Coefficient": corcoef,
        "R^2": corcoefSq,
        "Nash-Sutcliffe Coefficient": nsccoef,
        "lamda index": lamda_index,
    }

    return pd.Series(dic)
