import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np

def sgpvFun(results_df, null_interval, return_dg = True):
    """
    Calculates Second Generation P-values

    :param results_df: dataframe containing pyPheWAS results
    :param null_interval: null interval
    :param return_dg: return deltagap (gap between null interval & confidence interval normalized by null interval size)

    :type results_df: pandas DataFrame
    :type null_lb: numpy array (float)
    :type return_dg: boolean

    :returns: None
    """
    null_lb = null_interval[0]
    null_ub = null_interval[1]

    ci_size = results_df['ci.ub'] - results_df['ci.lb']
    null_size = null_ub - null_lb

    if (null_size == 0): # null interval is a single point
        null = null_lb
        f = lambda x: 0 if ((null < x['ci.lb']) | (null > x['ci.ub'])) else 0.5
        results_df['sgpv'] = results_df.apply(f,axis=1)

    else: # null interval is a range
        f = lambda x: max(null_lb,min(null_ub,x['ci.ub'])) - min(null_ub,max(null_lb,x['ci.lb']))
        overlap = results_df.apply(f,axis=1)
        results_df['sgpv'] = overlap/ci_size*np.maximum(ci_size/(2*null_size),1)

    if(return_dg):
        results_df['delta.gap'] = np.nan
        null = null_lb
        dg_exists = results_df['sgpv'] == 0
        if (null_size == 0): # null interval is a point
            f = lambda x: (x['ci.lb'] - null) if (null < x['ci.lb']) else (null - x['ci.ub'])
        else: # null interval is a range
            f = lambda x: (max(x['ci.lb'],null_lb) - min(null_ub,x['ci.ub']))/(null_size/2)
        results_df.loc[dg_exists,'delta.gap'] = results_df[dg_exists].apply(f,axis=1)

    return

def sgpvPowerFun(beta, beta_se, type = c("d", "c"), null_lb, null_ub):
    if(type == "d"):
        pnorm((null.lb)/beta.se - beta/beta.se - qnorm(1-0.05/2)) + pnorm(-(null.ub)/beta.se + beta/beta.se - qnorm(1-0.05/2))
    else if(type == "c"):
        pnorm((null.ub)/beta.se - beta/beta.se - qnorm(1-0.05/2)) - pnorm((null.lb)/beta.se - beta/beta.se + qnorm(1-0.05/2))

    return



def sgpvFdrFun(beta, beta_se, pi0 = 0.5, null, null_lb, null_ub):
    piA = 1 - pi0

    powerA = sgpvPowerFun(beta, beta.se, type = "d", null.lb = null.lb, null.ub = null.ub)
    power0 = sgpvPowerFun(null, beta.se, type = "d", null.lb = null.lb, null.ub = null.ub)

    fdr = (1 + powerA/power0 * piA/pi0)^(-1)
    return


def ppvFun(results_df, null_interval):
        """
        Calculates Positive Predictive Values

        :param results_df: dataframe containing pyPheWAS results
        :param null_interval: null interval

        :type results_df: pandas DataFrame
        :type null_lb: numpy array (float)

        :returns: None
        """
    ifelse(sgpv == 0, 1 - sgpvFdrFun(beta, beta.se, pi0 = pi0, null = point.null, null.lb(), null.ub()), NA))

    # if sgpv == 0,

    return
