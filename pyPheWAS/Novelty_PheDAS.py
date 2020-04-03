import numpy as np
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl


def secondGenPval(results_df, null_interval, return_dg=True):
    """
    Calculates Second Generation P-values

    :param results_df: dataframe containing pyPheWAS results
    :param null_interval: null interval
    :param return_dg: return deltagap (gap between null interval & confidence interval normalized by null interval size)

    :type results_df: pandas DataFrame
    :type null_interval: numpy array (float)
    :type return_dg: boolean

    :returns: None
    """
    null_lb = null_interval[0]
    null_ub = null_interval[1]

    ci_size = results_df['ci.ub'] - results_df['ci.lb']
    null_size = null_ub - null_lb

    if null_size == 0:  # null interval is a single point
        null = null_lb
        f = lambda x: 0 if ((null < x['ci.lb']) | (null > x['ci.ub'])) else 0.5
        results_df['sgpv'] = results_df.apply(f, axis=1)

    else:  # null interval is a range
        f = lambda x: max(null_lb, min(null_ub, x['ci.ub'])) - min(null_ub, max(null_lb, x['ci.lb']))
        overlap = results_df.apply(f, axis=1)
        results_df['sgpv'] = overlap / ci_size * np.maximum(ci_size / (2 * null_size), 1)

    if return_dg:
        results_df['delta.gap'] = np.nan
        null = null_lb
        dg_exists = results_df['sgpv'] == 0
        if null_size == 0:  # null interval is a point
            f = lambda x: (x['ci.lb'] - null) if (null < x['ci.lb']) else (null - x['ci.ub'])
        else:  # null interval is a range
            f = lambda x: (max(x['ci.lb'], null_lb) - min(null_ub, x['ci.ub'])) / (null_size / 2)
        results_df.loc[dg_exists, 'delta.gap'] = results_df[dg_exists].apply(f, axis=1)

    return


def sgpvPowerFun(beta, beta_se, null_lb, null_ub, ftype="d"):
    if ftype == "d":
        return norm.cdf(null_lb / beta_se - beta / beta_se - norm.ppf(1 - 0.05 / 2)) + \
               norm.cdf(-null_ub / beta_se + beta / beta_se - norm.ppf(1 - 0.05 / 2))
    elif ftype == "c":
        return norm.cdf(null_ub / beta_se - beta / beta_se - norm.ppf(1 - 0.05 / 2)) - \
               norm.cdf(null_lb / beta_se - beta / beta_se + norm.ppf(1 - 0.05 / 2))


def sgpvFdrFun(beta, beta_se, null, null_lb, null_ub, pi0=0.5):
    piA = 1 - pi0

    powerA = sgpvPowerFun(beta, beta_se, null_lb, null_ub, ftype="d")
    power0 = sgpvPowerFun(null, beta_se, null_lb, null_ub, ftype="d")

    fdr = (1 + powerA/power0 * piA/pi0) ** (-1)
    return fdr


def sgpvFcrFun(beta, beta_se, null, null_lb, null_ub, pi0=0.5):
    piA = 1 - pi0

    powerA = sgpvPowerFun(beta, beta_se, null_lb, null_ub, ftype="c")
    power0 = sgpvPowerFun(null, beta_se, null_lb, null_ub, ftype="c")

    fcr = (1 + power0/powerA * pi0/piA) ** (-1)
    return fcr


def positivePredictiveValue(results_df, null_interval, null_point, pi0=0.5):
    """
    Calculates Positive Predictive Value for each row in results_df with a second generation p-value equaling 0

    :param results_df: dataframe containing pyPheWAS results; required columns: beta, beta_se, sgpv
    :param null_interval: null interval
    :param null_point:
    :param pi0: null interval

    :type results_df: pandas DataFrame
    :type null_interval: numpy array (float)
    :type null_point: float
    :type pi0: float

    :returns: None
    """
    null_lb = null_interval[0]
    null_ub = null_interval[1]

    # if you test positive, how likely is it that you really have it?
    f = lambda x: (1 - sgpvFdrFun(x['beta'], x['beta.se'], null_point, null_lb, null_ub, pi0)) if x['sgpv'] == 0 else np.nan

    results_df['ppv'] = results_df.apply(f, axis=1)

    return


def noveltyScore(results_df):
    """
    Calculates Novelty Score for each row in results_df with a second generation p-value equaling 0

    :param results_df: dataframe containing pyPheWAS results; required columns: beta, beta_se, sgpv

    :type results_df: pandas DataFrame

    :returns: None
    """

    ecdfObs = ECDF(results_df['P.PM.phe.given.AD'])
    results_df['ecdf'] = ecdfObs(results_df['P.PM.phe.given.AD'])
    f = lambda x: (x['ppv'] * (1 - x['ecdf']) * 10) if x['sgpv'] == 0 else np.nan
    results_df['novelty_Findex'] = results_df.apply(f, axis=1)

    return

def plot_odds_ratio_novelty(regressions, null_interval, show_imbalance=True, save='', save_format=''):
    """
    Plots the data on a Log Odds Plot.

    :param regressions: dataframe containing the regression results
    :param save: the output file to save to (if empty, display the plot)
    :param show_imbalance: boolean variable that determines whether or not to show imbalances on the plot (default True)
    :param label_loc: the output file to save to (if empty, display the plot)
    :type regressions: pandas DataFrame
    :type save: str
    :type show_imbalance: boolean

    """

    # sort PheCodes by Novelty Index
    regressions.sort_values(by=["novelty_Findex"],inplace=True)

    # create colormap
    colormap = cm.get_cmap('plasma', 1000)

    # Initialize figure
    fig = plt.figure(1)
    ax = plt.subplot(111)
    frame1 = plt.gca()

    # Plot all points w/ labels
    e = 1  # vertical index
    text_size = 4
    axis_fontsize = 6
    artists = []
    phecode_labels = []
    phecode_locs = []
    plt.xlabel('Log odds ratio')
    for ix, data in regressions.iterrows():
        beta_ix = data['beta']
        # Add Phecode label
        if show_imbalance:
            if beta_ix > 0:
                artists.append(
                    ax.text(beta_ix+0.1, e+5, round(data['novelty_Findex'],3), rotation=0, ha='left', fontsize=text_size))
            else:
                artists.append(ax.text(beta_ix-0.1, e+5, round(data['novelty_Findex'],3), rotation=0, ha='right',
                                       fontsize=text_size))
        else:
            artists.append(
                ax.text(beta_ix+0.1, e+5, round(data['novelty_Findex'],3), rotation=0, va='bottom', fontsize=text_size))
        phecode_labels.append(data['PheName'])
        phecode_locs.append(e)

        # Plot Phecode Data
        ax.plot(beta_ix, e, 'o', fillstyle='full', markeredgewidth=0.0, color=colormap(data['novelty_Findex']/10.0))
        ax.plot([data['beta.lb'], data['beta.ub']], [e, e], color=colormap(data['novelty_Findex']/10.0))
        e += 15

    # Plot Null interval
    bottom, top = plt.ylim()
    rect = plt.Rectangle((null_interval[0], bottom), null_interval[1]-null_interval[0], top-bottom, facecolor="black", alpha=0.1)
    ax.add_patch(rect)

    # Plot y axis
    ax.axvline(x=0, color='black')
    plt.yticks(phecode_locs, phecode_labels, ha='right', fontsize=axis_fontsize)

    # Plot Colorbar
    # ax, _ = mpl.colorbar.make_axes(plt.gca(), shrink=0.5)
    # cbar = mpl.colorbar.ColorbarBase(ax, cmap=colormap, norm=mpl.colors.Normalize(vmin=0, vmax=10),)
    fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=10),cmap=colormap),
                 ax=ax,
                 ticks=range(0,11,2),
                 fraction=.1,
                 shrink=0.5,
                 label="Novelty Finding Index"
                 )


    # Save the plot
    if save:
        plt.savefig(save, format=save_format, bbox_extra_artists=artists, bbox_inches='tight', dpi=300)
        plt.clf()

    return