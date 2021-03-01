import numpy as np
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from ast import literal_eval
import pandas as pd
from tqdm import tqdm


def secondGenPval(results_df, null_interval, return_dg=False):
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

    ci_size = results_df['OR_uplim'] - results_df['OR_lowlim']
    null_size = null_ub - null_lb

    if null_size == 0:  # null interval is a single point
        null = null_lb
        f = lambda x: 0 if ((null < x['OR_lowlim']) | (null > x['OR_uplim'])) else 0.5
        results_df['sgpv'] = results_df.apply(f, axis=1)

    else:  # null interval is a range
        f = lambda x: max(null_lb, min(null_ub, x['OR_uplim'])) - min(null_ub, max(null_lb, x['OR_lowlim']))
        overlap = results_df.apply(f, axis=1)
        results_df['sgpv'] = overlap / ci_size * np.maximum(ci_size / (2 * null_size), 1)

    if return_dg:
        results_df['delta_gap'] = np.nan
        null = null_lb
        dg_exists = results_df['sgpv'] == 0
        if null_size == 0:  # null interval is a point
            f = lambda x: (x['OR_lowlim'] - null) if (null < x['OR_lowlim']) else (null - x['OR_uplim'])
        else:  # null interval is a range
            f = lambda x: (max(x['OR_lowlim'], null_lb) - min(null_ub, x['OR_uplim'])) / (null_size / 2)
        results_df.loc[dg_exists, 'delta_gap'] = results_df[dg_exists].apply(f, axis=1)

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
    f = lambda x: (1 - sgpvFdrFun(x['beta'], x['std_error'], null_point, null_lb, null_ub, pi0)) if x['sgpv'] == 0 else np.nan

    results_df['ppv'] = results_df.apply(f, axis=1)

    return


def calcNoveltyScore(results_df, null_interval):
    """
    Calculates Novelty Score for each row in results_df with a second generation p-value equaling 0

    :param results_df: dataframe containing pyPheWAS results

    :type results_df: pandas DataFrame

    :returns: None
    """
    # get the 2nd gen pvalues and positive predictive values
    secondGenPval(results_df, null_interval)
    positivePredictiveValue(results_df, null_interval, 1.0)

    # estimate the cumulative distribution function of the PubMed Proportions
    ecdfObs = ECDF(results_df['P_PM_phe_given_DX'])
    results_df['ecdf'] = ecdfObs(results_df['P_PM_phe_given_DX'])

    # combine everything to get the Novelty Finding Index
    f = lambda x: (x['ppv'] * (1 - x['ecdf']) * 10) if x['sgpv'] == 0 else np.nan
    results_df['Novelty_Finding_Index'] = results_df.apply(f, axis=1)

    return results_df


def get_joint_PubMed_articles(reg_df, dx_pubmed, pubmed_dir):
    print('Loading DX articles')
    # get count of DX articles on PubMed
    dx_set = set(literal_eval(dx_pubmed.loc[0, "IdsList"]))
    dx_count = len(dx_set)
    reg_df["DX_PM_count"] = dx_count

    # init other count columns
    reg_df["phe_PM_count"] = np.nan
    reg_df["joint_PM_count"] = np.nan
    reg_df["P_PM_phe_given_DX"] = np.nan

    print('Calculating PubMed Proportions')
    # iterate over all files output from the PubMed search
    file_list = [f for f in pubmed_dir.glob('phecode_pubmed_articles_*.csv')] # this list comprehension is unnecessary, but lets us easily count the number of files
    assert len(file_list) > 0, 'No PubMed search files found at %s' % pubmed_dir
    for pm_f in tqdm(file_list):
        phecode_pm = pd.read_csv(pubmed_dir / pm_f, dtype={"phewas_code": str})
        for ix, data in phecode_pm.iterrows():
            phe = data["phewas_code"] # change to PheWAS Code after testing
            if phe in reg_df["PheWAS Code"].values:
                phe_ix = reg_df["PheWAS Code"] == phe
                phe_set = set(literal_eval(data["IdsList"])) # get the set of article IDs found for this phecode
                joint_set = dx_set.intersection(phe_set) # get the set of article IDs that mention this phecode AND the DX
                reg_df.loc[phe_ix, "phe_PM_count"] = len(phe_set)
                reg_df.loc[phe_ix, "joint_PM_count"] = len(joint_set)
                if len(dx_set) != 0:
                    reg_df.loc[phe_ix, "P_PM_phe_given_DX"] = len(joint_set) / len(dx_set)

    return reg_df


def plot_log_odds_ratio_novelty(regressions, null_interval, save_f):
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
    regressions.sort_values(by=["Novelty_Finding_Index"],inplace=True)

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
    plt.xlabel('Log Odds Ratio')
    for ix, data in regressions.iterrows():
        log_odds = data['beta']
        # Add Novelty Finding Index annotation (only diff is text orientation & placement)
        if log_odds > 0:
            artists.append(ax.text(log_odds+0.1, e+5, round(data['Novelty_Finding_Index'],3),
                                   rotation=0, ha='left', fontsize=text_size))
        else:
            artists.append(ax.text(log_odds-0.1, e+5, round(data['Novelty_Finding_Index'],3),
                                   rotation=0, ha='right',fontsize=text_size))
        # Add Phenotype to label list for y-axis
        phecode_labels.append(data['PheWAS Name'])
        phecode_locs.append(e)
        # Plot Phecode Data
        ax.plot(log_odds, e, 'o', fillstyle='full', markeredgewidth=0.0, color=colormap(data['Novelty_Finding_Index']/10.0))
        ax.plot([data['beta_lowlim'], data['beta_uplim']], [e, e], color=colormap(data['Novelty_Finding_Index']/10.0))
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
                 ax=ax, ticks=range(0,11,2), fraction=.1,
                 shrink=0.5, label="Novelty Finding Index"
                 )
    # Save the plot
    plt.savefig(save_f, bbox_extra_artists=artists, bbox_inches='tight', dpi=300)
    plt.clf()

    return