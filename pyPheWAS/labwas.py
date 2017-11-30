output_columns = ['PheWAS Code',
                  'PheWAS Name',
                  'p-val',
                  '\"-log(p)\"',
                  'beta',
                  'Conf-interval beta',
                  'cpt']

imbalance_colors = {
    0: 'white',
    1: 'deepskyblue',
    -1: 'red'
}
m = len(fm[0])
p_values = np.zeros(m, dtype=float)
icodes = []
# store all of the pertinent data from the regressions
regressions = pd.DataFrame(columns=output_columns)
labnames=df.columns

def get_bon_thresh(normalized, power):  # same
    """
    Calculate the bonferroni correction threshold.

    Divide the power by the sum of all finite values (all non-nan values).

    :param normalized: an array of all normalized p-values. Normalized p-values are -log10(p) where p is the p-value.
    :param power: the threshold power being used (usually 0.05)
    :type normalized: numpy array
    :type power: float

    :returns: The bonferroni correction
    :rtype: float

    """
    return power / sum(np.isfinite(normalized))

for index in range(m):
    print index
    phen_vector1 = fm[:, index]
    res = calculate_odds_ratio(genotypes, phen_vector1,0)

    # save all of the regression data
    phewas_info = [labnames[index],labnames[index],labnames[index]]
    stat_info = res[2]
    info = phewas_info[0:2] + [res[1]] + stat_info + [phewas_info[2]]
    regressions.loc[index] = info

    p_values[index] = res[1]


def get_imbalances(regressions):
    """
    Generates a numpy array of the imbalances.

    For a value *x* where *x* is the beta of a regression:

    ========= ====== =======================================================
    *x* < 0   **-1** The regression had a negative beta value
    *x* = nan **0**  The regression had a nan beta value (and a nan p-value)
    *x* > 0   **+1** The regression had a positive beta value
    ========= ====== =======================================================

    These values are then used to get the correct colors using the imbalance_colors.

    :param regressions: DataFrame containing a variety of different output values from the regression performed. The only one used for this function are the 'beta' values.
    :type regressions: pandas DataFrame

    :returns: A list that is the length of the number of regressions performed. Each element in the list is either a -1, 0, or +1. These are used as explained above.
    :rtype: numpy array
    """

    imbalance = np.array(regressions['beta'])
    imbalance[np.isnan(imbalance)] = 0
    imbalance[imbalance > 0] = 1
    imbalance[imbalance < 0] = -1
    return imbalance

def calculate_odds_ratio(genotypes, phen_vector1,reg_type):  # diff - done

    data = genotypes
    data['y'] = phen_vector1
    f = 'genotype ~ y'
    try:
        if reg_type == 0:
            logreg = smf.logit(f, data).fit(method='bfgs', disp=False)
            p = logreg.pvalues.y
            odds = logreg.params.y
            conf = logreg.conf_int()
            od = [-math.log10(p), logreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
        else:
            linreg = smf.logit(f, data).fit(method='bfgs', disp=False)
            p = linreg.pvalues.y
            odds = linreg.params.y
            conf = linreg.conf_int()
            od = [-math.log10(p), linreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
    except:
        odds = 0
        p = np.nan
        od = [np.nan, np.nan, np.nan]
    return (odds, p, od)

def plot_data_points(y, thresh,labnames,save='', imbalances=np.array([])):  # same

    idx = y.sort_values().index

    # Plot each of the points, if necessary, label the points.
    e = 1
    artists = []
    for i in idx:
        if imbalances[i] >0:
            plt.plot(e, y[i], 'o', color=imbalance_colors[imbalances[i]], fillstyle='full', markeredgewidth=0.0)
        if y[i] > thresh and imbalances[i] > 0:
            artists.append(plt.text(e, y[i], labnames[i], fontsize=5,rotation=70, va='bottom'))
            e += 10

    # If the imbalance is to be shown, draw lines to show the categories.
    # if show_imbalance:
    #     for pos in linepos:
    #         plt.axvline(x=pos, color='black', ls='dotted')

    # Plot a blue line at p=0.05 and plot a red line at the line for the threshold type.
    plt.axhline(y=-math.log10(0.05), color='blue')
    plt.axhline(y=thresh, color='red')

    # Set windows and labels
    # plt.xticks(x_label_positions, x_labels, rotation=70, fontsize=10)
    plt.ylim(ymin=0, ymax=max(y[imbalances>0])+5)
    plt.xlim(xmin=0, xmax=e)
    plt.ylabel('-log10(p)')

    # Determine the type of output desired (saved to a plot or displayed on the screen)
    if save:
        pdf = PdfPages(save)
        pdf.savefig(bbox_extra_artists=artists, bbox_inches='tight')
        pdf.close()
    else:
        plt.subplots_adjust(left=0.05, right=0.85)
        plt.show()

    # Clear the plot in case another plot is to be made.
    plt.clf()


regressions[(y > -math.log10(0.05))&(imbalances<0)].to_csv('labwasneg.csv')