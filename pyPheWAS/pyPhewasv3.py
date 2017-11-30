# Shikha Chaganti
# Kunal Nabar
# Vanderbilt University
# Medical-image Analysis and Statistical Interpretation Lab
# newphewas
# v2.0

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter
import time, math, scipy.stats
import pandas as pd
from matplotlib import rcParams
import os
import statsmodels.formula.api as smf
import matplotlib.lines as mlines
import statsmodels.discrete.discrete_model as sm

"""
I/O Reading Input From Files
"""


def get_codes():  # same
    """
	Gets the PheWAS codes from a local csv file and load it into a pandas DataFrame.

	:returns: All of the codes from the resource file.
	:rtype: pandas DataFrame

	"""
    sep = os.sep
    path = os.path.dirname(os.path.abspath(__file__))
    filename = os.sep.join([path, 'resources', 'codes.csv'])
    return pd.read_csv(filename)


def get_input(path, filename,reg_type):  # diff -done - add duration
    """
	Read all of the phenotype data from the given file and load it into a pandas DataFrame.

	:param path: The path to the file that contains the phenotype data
	:param filename: The name of the file that contains the phenotype data.
	:type path: string
	:type filename: string

	:returns: The data from the phenotype file.
	:rtype: pandas DataFrame
	"""
    wholefname = path + filename
    icdfile = pd.read_csv(wholefname)
    icdfile['icd9'] = icdfile['icd9'].str.strip()
    # if reg_type == 0:
    #     # g=icdfile.groupby(['id','icd9'])
    #     phenotypes = pd.merge(icdfile, codes, on='icd9')
    #     phenotypes['MaxAgeAtICD'] = 0
    #     phenotypes['MaxAgeAtICD'] = phenotypes.groupby(['id', 'phewas_code'])['AgeAtICD9'].transform('max')
    # else:
    """
    This needs to be changed, need to adjust for a variety of different naming conventions
    in the phenotype file, not simply 'AgeAtICD', 'id', 'icd9', etc.
    Either we need to adjust for different names in the code, or state explicitly in the
    documentation that we cannot do things like this.
    """
    phenotypes = pd.merge(icdfile, codes, on='icd9')
    phenotypes['count'] = 0
    phenotypes['count'] = phenotypes.groupby(['id', 'phewas_code'])['count'].transform('count')
    phenotypes['MaxAgeAtICD'] = 0
    phenotypes['MaxAgeAtICD'] = phenotypes.groupby(['id', 'phewas_code'])['AgeAtICD9'].transform('max')
    phenotypes['duration'] = phenotypes.groupby(['id', 'phewas_code'])['AgeAtICD9'].transform('max') - \
                             phenotypes.groupby(['id', 'phewas_code'])['AgeAtICD9'].transform('min') + 1
    phenotypes['lor'] = phenotypes.groupby('id')['AgeAtICD9'].transform('max') - \
                        phenotypes.groupby('id')['AgeAtICD9'].transform('min') + 1
    return phenotypes


def get_phewas_info(p_index):  # same
    """
	Returns all of the info of the phewas code at the given index.

	:param p_index: The index of the desired phewas code
	:type p_index: int

	:returns: A list including the code, the name, and the rollup of the phewas code. The rollup is a list of all of the ICD-9 codes that are grouped into this phewas code.
	:rtype: list of strings
	"""
    p_code = phewas_codes.loc[p_index].phewas_code
    corresponding = codes[codes.phewas_code == p_code]

    p_name = corresponding.iloc[0].phewas_string
    p_rollup = ','.join(codes[codes.phewas_code == p_code].icd9.tolist())
    return [p_code, p_name, p_rollup]


def get_group_file(path, filename):  # same
    """
	Read all of the genotype data from the given file and load it into a pandas DataFrame.

	:param path: The path to the file that contains the phenotype data
	:param filename: The name of the file that contains the phenotype data.
	:type path: string
	:type filename: string

	:returns: The data from the genotype file.
	:rtype: pandas DataFrame
	"""
    wholefname = path + filename
    genotypes = pd.read_csv(wholefname)
    return genotypes


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


def generate_feature_matrix(genotypes, phenotypes, reg_type,phewas_cov=''):  # diff - done
    """
	Generates the feature matrix that will be used to run the regressions.

	:param genotypes:
	:param phenotypes:
	:type genotypes:
	:type phenotypes:

	:returns:
	:rtype:

	"""
    feature_matrix = np.zeros((3, genotypes.shape[0], phewas_codes.shape[0]), dtype=float)
    count = 0
    for i in genotypes['id']:
        if reg_type == 0:
            temp = pd.DataFrame(phenotypes[phenotypes['id'] == i][['phewas_code', 'MaxAgeAtICD','count']]).drop_duplicates()
            match = phewas_codes['phewas_code'].isin(list(phenotypes[phenotypes['id'] == i]['phewas_code']))
            cts = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['count']
            cts[np.isnan(cts)] = 0
            match = (match)&(cts>0)
            feature_matrix[0][count, match[match == True].index] = 1

            age = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['MaxAgeAtICD']
            #assert np.all(np.isfinite(age)), "make sure MaxAgeAtVisit is filled"
            age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeBeforeDx']
            feature_matrix[1][count, :] = age
            if phewas_cov:
                feature_matrix[2][count, :] = int(phewas_cov in list(phenotypes[phenotypes['id'] == i]['phewas_code']))


        else:
            if reg_type == 1:
                temp = pd.DataFrame(
                    phenotypes[phenotypes['id'] == i][['phewas_code', 'MaxAgeAtICD', 'count','lor']]).drop_duplicates()
                cts = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['count']
                cts[np.isnan(cts)] = 0
                if temp.empty!=1:
                    cts=cts/temp['lor'].iloc[0]
                feature_matrix[0][count, :] = cts
                age = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['MaxAgeAtICD']
                #assert np.all(np.isfinite(age)), "make sure MaxAgeAtVisit is filled"
                age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeBeforeDx']
                feature_matrix[1][count, :] = age
                if phewas_cov:
                    feature_matrix[2][count, :] = int(
                        phewas_cov in list(phenotypes[phenotypes['id'] == i]['phewas_code']))
            elif reg_type == 2:
                temp = pd.DataFrame(
                    phenotypes[phenotypes['id'] == i][['phewas_code', 'MaxAgeAtICD', 'duration','lor']]).drop_duplicates()
                dura = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['duration']
                dura[np.isnan(dura)] = 0
                if temp.empty!=1:
                    dura=dura/temp['lor'].iloc[0]
                feature_matrix[0][count, :] = dura
                age = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['MaxAgeAtICD']
                #assert np.all(np.isfinite(age)), "make sure MaxAgeAtVisit is filled"
                age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeBeforeDx']
                feature_matrix[1][count, :] = age
                if phewas_cov:
                    feature_matrix[2][count, :] = int(
                        phewas_cov in list(phenotypes[phenotypes['id'] == i]['phewas_code']))

        count += 1
    return feature_matrix


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


def get_fdr_thresh(p_values, power):
    """
	Calculate the false discovery rate threshold.

	:param p_values: a list of p-values obtained by executing the regression
	:param power: the thershold power being used (usually 0.05)
	:type p_values: numpy array
	:type power: float

	:returns: the false discovery rate
	:rtype: float
	"""
    sn = np.sort(p_values)
    sn = sn[np.isfinite(sn)]
    sn = sn[::-1]
    for i in range(len(sn)):
        thresh = power * i / len(sn)
        if sn[i] <= thresh:
            break
    return sn[i]

def get_bhy_thresh(p_values, power):
    """
	Calculate the false discovery rate threshold.

	:param p_values: a list of p-values obtained by executing the regression
	:param power: the thershold power being used (usually 0.05)
	:type p_values: numpy array
	:type power: float

	:returns: the false discovery rate
	:rtype: float
	"""
    sn = np.sort(p_values)
    sn = sn[np.isfinite(sn)]
    sn = sn[::-1]
    for i in range(len(sn)):
        thresh = power * i / (8.1*len(sn))
        if sn[i] <= thresh:
            break
    return sn[i]



def run_phewas(fm, genotypes, covariates, reg_type, response='', phewas_cov=''):  # same
    """
	For each phewas code in the feature matrix, run the specified type of regression and save all of the resulting p-values.

	:param fm: The phewas feature matrix.
	:param genotypes: A pandas DataFrame of the genotype file.
	:param covariates: The covariates that the function is to be run on.

	:returns: A tuple containing indices, p-values, and all the regression data.
	"""
    m = len(fm[0, 0])
    p_values = np.zeros(m, dtype=float)
    icodes = []
    # store all of the pertinent data from the regressions
    regressions = pd.DataFrame(columns=output_columns)
    control = fm[0][genotypes.genotype == 0, :]
    disease = fm[0][genotypes.genotype == 1, :]
    inds = np.where((control.any(axis=0) & ~disease.any(axis=0)) | (~control.any(axis=0) & disease.any(axis=0)))[0]
    for index in range(m):
        phen_vector1 = fm[0][:, index]
        phen_vector2 = fm[1][:, index]
        phen_vector3 = fm[2][:, index]
        if sum(phen_vector1)>5:
            if index in inds:
                print index
                res = calculate_odds_ratio(genotypes, phen_vector1, phen_vector2, reg_type, covariates, lr=1, response=response,
                                   phen_vector3=phen_vector3)
            else:
                res = calculate_odds_ratio(genotypes, phen_vector1, phen_vector2, reg_type, covariates, lr=0,
                                           response=response,
                                           phen_vector3=phen_vector3)
        else:
            odds = 0
            p = 1
            od = [-0.0, 1.0, 0.0, np.nan]
            res = (odds, p, od)

        # save all of the regression data

        phewas_info = get_phewas_info(index)
        stat_info = res[2]
        info = phewas_info[0:2] +stat_info + [phewas_info[2]]
        regressions.loc[index] = info

        p_values[index] = res[1]
    return (np.array(range(m)), p_values, regressions)


"""
Plotting
"""


def get_x_label_positions(categories, lines=True):  # same
    """
	This method is used get the position of the x-labels and the lines between the columns

	:param categories: list of the categories
	:param lines: a boolean which determines the locations returned (either the center of each category or the end)
	:type categories:
	:type lines: bool

	:returns: A list of positions
	:rtype: list of ints

	"""
    tt = Counter(categories)
    s = 0
    label_positions = []
    for _, v in tt.items():
        if lines:
            inc = v // 2
        else:
            inc = v
        label_positions.append(s + inc)
        s += v
    return label_positions


def plot_data_points(x, y, thresh0, thresh1, thresh2, thresh_type, save='', path='', imbalances=np.array([])):  # same
    """
    Plots the data with a variety of different options.

    This function is the primary plotting function for pyPhewas.

    :param x: an array of indices
    :param y: an array of p-values
    :param thresh: the threshold power
    :param save: the output file to save to (if empty, display the plot)
    :param imbalances: a list of imbalances
    :type x: numpy array
    :type y: numpy array
    :type thresh: float
    :type save: str
    :type imbalances: numpy array

    """

    # Determine whether or not to show the imbalance.
    fig = plt.figure()
    ax = plt.subplot(111)
    show_imbalance = imbalances.size != 0

    # Sort the phewas codes by category.
    c = codes.loc[phewas_codes['index']]
    c = c.reset_index()
    idx = c.sort_values(by='category').index

    # Get the position of the lines and of the labels
    # linepos = get_x_label_positions(c['category'].tolist(), False)
    # x_label_positions = get_x_label_positions(c['category'].tolist(), True)
    # x_labels = c.sort_values('category').category_string.drop_duplicates().tolist()

    # Plot each of the points, if necessary, label the points.
    e = 1
    artists = []
    frame1 = plt.gca()
    # ax.axhline(y=-math.log10(0.05), color='yellow', ls='dotted')
    ax.axhline(y=thresh0, color='red', ls='dotted')
    ax.axhline(y=thresh1, color='yellow', ls='dotted')
    ax.axhline(y=thresh2, color='orange', ls='dotted')
    # ax.xticks(x_label_positions, x_labels, rotation=70, fontsize=10)
    # ax.xlim(xmin=0, xmax=len(c))
    plt.ylabel('-log10(p)')
    if thresh_type == 0:
        thresh = thresh0
    elif thresh_type == 1:
        thresh = thresh1
    else:
        thresh = thresh2

    y_label_positions = [thresh0, thresh1,thresh2]
    plt.yticks(y_label_positions, ['Bonf p = ' + '{:.2e}'.format(np.power(10, -thresh0)),
                                   'Benj-Hoch p = ' + str(round(np.power(10, -thresh1), 3)),
                                'Benj-Hoch-Yek p = ' + str(round(np.power(10, -thresh2), 3))], rotation=10, fontsize=10)

    for i in idx:
        if y[i] > thresh:
            e += 15
            if show_imbalance:  # and imbalances[i]>0:
                # if imbalances[i]>0:
                artists.append(ax.text(e, y[i], c['phewas_string'][i], rotation=89, va='bottom', fontsize=8))
            # else:
            #	artists.append(ax.text(e, -y[i], c['phewas_string'][i], rotation=271, va='top',fontsize=8))
            elif not show_imbalance:
                artists.append(ax.text(e, y[i], c['phewas_string'][i], rotation=40, va='bottom'))
        else:
            e += 0

        if show_imbalance:
            if y[i] > thresh:
                if imbalances[i] > 0:
                    ax.plot(e, y[i], '+', color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full',
                            markeredgewidth=1.5)

                else:
                    # ax.plot(e,y[i],'o', color=plot_colors[c[i:i+1].category_string.values[0]], fillstyle='full', markeredgewidth=0.0)
                    ax.plot(e, y[i], '_', color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full',
                            markeredgewidth=1.5)
        else:
            ax.plot(e, y[i], 'o', color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full',
                    markeredgewidth=0.0)
    line1 = []
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.05, box.width, box.height*0.95])
    for lab in plot_colors.keys():
        line1.append(
            mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=plot_colors[lab], label=lab))
    artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, 0), loc='upper center', fancybox=True, ncol=4, prop={'size': 6}))
    ax.axhline(y=0, color='black')
    frame1.axes.get_xaxis().set_visible(False)

    # If the imbalance is to be shown, draw lines to show the categories.
    # if show_imbalance:
    # 	for pos in linepos:
    # 		ax.axvline(x=pos, color='black', ls='dotted')
    # Determine the type of output desired (saved to a plot or displayed on the screen)
    if save:
        pdf = PdfPages(path + save)
        pdf.savefig(bbox_extra_artists=artists, bbox_inches='tight')
        pdf.close()
    else:
        ax.subplots_adjust(left=0.05, right=0.85)
        ax.show()

    # Clear the plot in case another plot is to be made.
    plt.clf()


def plot_odds_ratio(y, p, thresh0, thresh1, thresh2, thresh_type, save='', path='', imbalances=np.array([])):  # same
    """
	Plots the data with a variety of different options.

	This function is the primary plotting function for pyPhewas.

	:param x: an array of indices
	:param y: an array of p-values
	:param thresh: the threshold power
	:param save: the output file to save to (if empty, display the plot)
	:param imbalances: a list of imbalances
	:type x: numpy array
	:type y: numpy array
	:type thresh: float
	:type save: str
	:type imbalances: numpy array

	"""

    # Determine whether or not to show the imbalance.
    fig = plt.figure()
    ax = plt.subplot(111)
    show_imbalance = imbalances.size != 0

    # Sort the phewas codes by category.
    c = codes.loc[phewas_codes['index']]
    c = c.reset_index()
    idx = c.sort_values(by='category').index

    # Get the position of the lines and of the labels
    # linepos = get_x_label_positions(c['category'].tolist(), False)
    # x_label_positions = get_x_label_positions(c['category'].tolist(), True)
    # x_labels = c.sort_values('category').category_string.drop_duplicates().tolist()

    # Plot each of the points, if necessary, label the points.
    e = 1
    artists = []
    frame1 = plt.gca()
    # ax.xticks(x_label_positions, x_labels, rotation=70, fontsize=10)
    plt.xlabel('Log odds ratio')

    if thresh_type == 0:
        thresh = thresh0
    elif thresh_type == 1:
        thresh = thresh1
    else:
        thresh = thresh2

    # plt.xlim(xmin=min(y[p>thresh,1]), xmax=max(y[p>thresh,2]))

    for i in idx:
        if p[i] > thresh:
            e += 15
            if show_imbalance:  # and imbalances[i]>0:
                if imbalances[i] > 0:
                    artists.append(ax.text(y[i][0], e, c['phewas_string'][i], rotation=0, ha='left', fontsize=6))
                else:
                    artists.append(ax.text(y[i][0], e, c['phewas_string'][i], rotation=0, ha='right', fontsize=6))
            elif not show_imbalance:
                artists.append(ax.text(e, y[i][0], c['phewas_string'][i], rotation=40, va='bottom'))
        else:
            e += 0

        if show_imbalance:
            if p[i] > thresh:
                ax.plot(y[i][0], e, 'o', color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full',
                        markeredgewidth=0.0)
                ax.plot([y[i, 1], y[i, 2]], [e, e], color=plot_colors[c[i:i + 1].category_string.values[0]])
            # else:
            # ax.plot(e,y[i],'o', color=plot_colors[c[i:i+1].category_string.values[0]], fillstyle='full', markeredgewidth=0.0)
            #	ax.plot(e,-y[i],'o', color=plot_colors[c[i:i+1].category_string.values[0]], fillstyle='full', markeredgewidth=0.0)
        else:
            ax.plot(e, y[i], 'o', color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full',
                    markeredgewidth=0.0)
    line1 = []
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height*0.05, box.width, box.height*0.95])
    for lab in plot_colors.keys():
        line1.append(
            mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=plot_colors[lab], label=lab))
    artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, -0.15), loc='upper center', fancybox=True, ncol=4, prop={'size': 6}))
    ax.axvline(x=0, color='black')
    frame1.axes.get_yaxis().set_visible(False)

    # If the imbalance is to be shown, draw lines to show the categories.
    # if show_imbalance:
    # 	for pos in linepos:
    # 		ax.axvline(x=pos, color='black', ls='dotted')
    # Determine the type of output desired (saved to a plot or displayed on the screen)
    if save:
        pdf = PdfPages(path + save)
        pdf.savefig(bbox_extra_artists=artists, bbox_inches='tight')
        pdf.close()
    else:
        ax.subplots_adjust(left=0.05, right=0.85)
        ax.show()

    # Clear the plot in case another plot is to be made.
    plt.clf()


def calculate_odds_ratio(genotypes, phen_vector1, phen_vector2, reg_type, covariates, lr=0, response='',
                         phen_vector3=''):  # diff - done
    """
	Runs the regression for a specific phenotype vector relative to the genotype data and covariates.

	:param genotypes: a DataFrame containing the genotype information
	:param phen_vector: a array containing the phenotype vecto
	:param covariates: a string containing all desired covariates
	:type genotypes: pandas DataFrame
	:type phen_vector: numpy array
	:type covariates: string

	.. note::
		The covariates must be a string that is delimited by '+', not a list.
		If you are using a list of covariates and would like to convert it to the pyPhewas format, use the following::

			l = ['genotype', 'age'] # a list of your covariates
			covariates = '+'.join(l) # pyPhewas format

		The covariates that are listed here *must* be headers to your genotype CSV file.
	"""

    data = genotypes
    data['y'] = phen_vector1
    data['MaxAgeAtICD'] = phen_vector2

    # data[data['SEX']>0]=0
    if response:
        f = response + '~ genotype + ' + covariates
        if phen_vector3.any():
            data['phe'] = phen_vector3
            f = response + '~ genotype + phe +' + covariates
    else:
        f = 'genotype~y+' + covariates
        if phen_vector3.any():
            data['phe'] = phen_vector3
            f = 'genotype ~ phe +' + covariates
    try:
        #if reg_type == 0:
        if lr >2:
            logreg = smf.logit(f, data).fit(disp=False)
            # logit = sm.Logit(data['genotype'], data[['y', 'MaxAgeAtVisit', 'sex']])
            # lf = logit.fit_regularized(method='l1', alpha=0.9,disp=0,trim_mode='size',qc_verbose=0)
            p = logreg.pvalues.y
            # p = lf.pvalues.y
            odds = 0  # logreg.deviance
            conf = logreg.conf_int()
            # conf = lf.conf_int()
            od = [-math.log10(p), p, logreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
            # od = [-math.log10(p), p, lf.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
        # od=[np.nan,np.nan,np.nan]
        #elif reg_type > 3:
        elif lr == 10:
            # linreg = smf.logit(f, data).fit(disp=False)
            logit = sm.Logit(data['genotype'], data[['y', 'MaxAgeAtVisit', 'sex']])
            lf = logit.fit_regularized(method='l1', alpha=1, disp=0,trim_mode='size',qc_verbose=0)
            # p = linreg.pvalues.y
            p = lf.pvalues.y
            odds = 0
            # conf = linreg.conf_int()
            conf = lf.conf_int()
            # od = [-math.log10(p), p, linreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
            od = [-math.log10(p), p, lf.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
        else:
            linreg = smf.logit(f, data).fit(method='bfgs', disp=False)
            # logit = sm.Logit(data['genotype'], data[['y', 'MaxAgeAtVisit', 'sex']])
            # lf = logit.fit_regularized(method='l1', alpha=0.7, disp=0,trim_mode='size',qc_verbose=0)
            p = linreg.pvalues.y
            # p = lf.pvalues.y
            odds = 0
            conf = linreg.conf_int()
            # conf = lf.conf_int()
            od = [-math.log10(p), p, linreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
            # od = [-math.log10(p), p, lf.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
    except:
        odds = 0
        p = np.nan
        od = [np.nan, np.nan, np.nan, np.nan]
    return (odds, p, od)


"""
Begin init code
"""
test = 1
codes = get_codes()
phewas_codes = pd.DataFrame(codes['phewas_code'].drop_duplicates());
phewas_codes = phewas_codes.reset_index()

output_columns = ['PheWAS Code',
                  'PheWAS Name',
                  '\"-log(p)\"',
                  'p-val',
                  'beta',
                  'Conf-interval beta',
                  'ICD-9']
plot_colors = {'-': 'gold',
               'circulatory system': 'red',
               'congenital anomalies': 'mediumspringgreen',
               'dermatologic': 'seagreen',
               'digestive': 'yellowgreen',
               'endocrine/metabolic': 'darkred',
               'genitourinary': 'darkkhaki',
               'hematopoietic': 'orange',
               'infectious diseases': 'blue',
               'injuries & poisonings': 'slategray',
               'mental disorders': 'fuchsia',
               'musculoskeletal': 'darkgreen',
               'neoplasms': 'teal',
               'neurological': 'olive',
               'pregnancy complications': 'peachpuff',
               'respiratory': 'brown',
               'sense organs': 'darkviolet',
               'symptoms': 'aqua'}
imbalance_colors = {
    0: 'white',
    1: 'deepskyblue',
    -1: 'red'
}

gen_ftype = 0
neglogp = np.vectorize(lambda x: -math.log10(x) if x != 0 else 0)


def phewas(path, filename, groupfile, covariates, response='', phewas_cov='', reg_type=0, thresh_type=0, control_age=0,
           save='', saveb='', output='', show_imbalance=False):  # same
    """
	The main phewas method. Takes a path, filename, groupfile, and a variety of different options.

	:param path: the path to the file that contains the phenotype data
	:param filename: the name of the phenotype file.
	:param groupfile: the name of the genotype file.
	:param covariates: a list of covariates.
	:param reg_type: the type of regression to be used
	:param thresh_type: the type of threshold to be used
	:param save: the desired filename to save the phewas plot
	:param output: the desired filename to save the regression output
	:param show_imbalance: determines whether or not to show the imbalance
	:type path: st
	:type filename: str
	:type groupfile: str
	:type covariates: str
	:type reg_type: int
	:type thresh_type: int
	:type save: str
	:type output: str
	:type show_imbalance: bool
	"""
    start_time = time.time()
    global codes, phewas_codes, gen_ftype, neglogp

    print("reading in data")

    gen_ftype = reg_type
    phenotypes = get_input(path, filename,reg_type)
    genotypes = get_group_file(path, groupfile)
    fm = generate_feature_matrix(genotypes, phenotypes, reg_type, phewas_cov)
    # print(len(fm))
    if response:
        results = run_phewas(fm, genotypes, covariates, reg_type, response=response, phewas_cov=phewas_cov)
    else:
        results = run_phewas(fm, genotypes, covariates, reg_type, phewas_cov=phewas_cov)

    regressions = results[2]

    normalized = neglogp(results[1])
    # if thresh_type==0:
    thresh0 = get_bon_thresh(normalized, 0.05)
    # elif thresh_type==1:
    thresh1 = get_fdr_thresh(results[1], 0.05)
    thresh2 = get_bhy_thresh(results[1], 0.05)
    imbalances = np.array([])
    if show_imbalance:
        imbalances = get_imbalances(regressions)
    plot_data_points(results[0], normalized, -math.log10(thresh0), -math.log10(thresh1),-math.log10(thresh2), thresh_type, save, path,
                     imbalances)

    regressions[['lowlim', 'uplim']] = regressions['Conf-interval beta'].str.split(',', expand=True)
    regressions.uplim = regressions.uplim.str.replace(']', '')
    regressions.lowlim = regressions.lowlim.str.replace('[', '')
    y = regressions[['beta', 'lowlim', 'uplim']].as_matrix()
    y = y.astype(float)
    plot_odds_ratio(y, normalized, -math.log10(thresh0),  -math.log10(thresh1),-math.log10(thresh2), thresh_type, saveb, path, imbalances)

    sig_regressions = regressions.dropna(subset=['"-log(p)"']).sort_values('"-log(p)"', ascending=False)
    if thresh_type == 0:
        sig_regressions = sig_regressions[sig_regressions['"-log(p)"'] > -math.log10(thresh0)]
    else:
        sig_regressions = sig_regressions[sig_regressions['"-log(p)"'] > -math.log10(thresh1)]
    if output:
        sig_regressions.to_csv(path + output, index=False)
    return (results[0], results[1], regressions)

