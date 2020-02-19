#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import numpy as np


def hind_mse(obsrv,model):
#    Mean Absolute Error
   
    if len(obsrv)==len(model):
        mse = np.nanmean(abs(obsrv-model))
    else:
        print('dimensions mismatch')
    return mse



def hind_rmse(obsrv,model):
#    RMSE
   
    if len(obsrv)==len(model):
        rmse = np.sqrt(np.nanmean((obsrv-model)**2))
    else:
        print('dimensions mismatch')
    return rmse


def hind_scatterIndex(obsrv,model):
#    Scatter Index 
   
    if len(obsrv)==len(model):
        scInd = np.sqrt(np.nanmean((obsrv-model)**2))/np.nanmean(obsrv)
    else:
        print('dimensions mismatch')
    return scInd



def hind_perc_rmse(obsrv,model):
#   percentage RMSE

    if len(obsrv)==len(model):
        percrmse = 100*(hind_rmse(obsrv,model))/(np.nanmax(obsrv))
    else:
        print('dimensions mismatch')
    return percrmse


def hind_bias(obsrv,model):
#   BIAS or mean error
# positive values indicate positively biased modelled data
# negative values indicate negatively biased modelled data

    if len(obsrv)==len(model):
        bias=np.nanmean(model-obsrv)
    else:
        print('dimensions mismatch')
    return bias


def hind_sdr(obsrv,model):
# Standard deviation of residuals
# reference: https://cirpwiki.info/wiki/Statistics
    
    if len(obsrv)==len(model):
        sdr = np.sqrt(np.nanmean((model-obsrv-np.nanmean(model)+np.nanmean(obsrv))**2))
    else:
        print('dimensions mismatch')
    return sdr


def hind_corcoef(obsrv,model):
# Standard deviation of residuals
# reference: https://cirpwiki.info/wiki/Statistics
    s1N=obsrv[~np.isnan(obsrv)]
    s2N=model[~np.isnan(model)]
    
    if len(s1N)==len(s2N):
        corcoef = np.corrcoef(s1N,s2N)[0,1]
    else:
        print('dimensions mismatch')
    return corcoef


def hind_Rsquared(obsrv,model):
# R**2
# reference: https://cirpwiki.info/wiki/Statistics
    s1N=s1[~np.isnan(obsrv)]
    s2N=s1[~np.isnan(model)]
    
    if len(S1N)==len(s2N):
        corcoefSq = np.corrcoef(s1N,s2N)[0,1]**2
    else:
        print('dimensions mismatch')
    return corcoefSq


def hind_Nash_Sutcliffe(obsrv,model):
# Nash-Sutcliffe Coefficient (Nash and Sutcliffe 1970)
#Range       Qualification
#0.8<PS<1.0  Excellent
#0.6<PS<0.8  Good
#0.3<PS<0.6  Reasonable
#0<PS<0.3    Poor
#PS<0        Bad 
# reference: https://cirpwiki.info/wiki/Statistics

    if len(obsrv)==len(model):
        nsccoef = 1 - np.nanmean((model-obsrv)**2)/np.nanmean((obsrv-np.nanmean(obsrv))**2)
    else:
        print('dimensions mismatch')
    return nsccoef

 
    
def hind_lamda(obsrv,model):
# lamda index of Duveiller et al (2016) in ScRe
# https://www.nature.com/articles/srep19401
# example of application in R
# http://eoscience.esa.int/landtraining2017/files/materials/D5P3a_I.pdf
    
    if len(obsrv)==len(model) & len(obsrv)>1:
        Xmean=np.nanmean(obsrv)
        Ymean=np.nanmean(model)
        nObs=len(obsrv)
        if hind_corcoef(obsrv,model)>=0:
            kappa=0
            #print('corcoef is positive',kappa)
        else:
            kappa=2*abs(np.nansum((obsrv-Xmean)*(model-Ymean)))
            #print('corcoef is negative',kappa)
        
        Nomin=np.nansum((obsrv-model)**2)
        Denom=np.nansum((obsrv-Xmean)**2)+np.nansum((model-Ymean)**2)+nObs*((Xmean-Ymean)**2)+kappa
        lamda_index=1-Nomin/Denom
    else:
        print('dimensions mismatch')
    return lamda_index
        










