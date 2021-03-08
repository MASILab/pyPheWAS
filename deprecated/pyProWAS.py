from collections import Counter
import getopt
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import pandas as pd
import scipy.stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from tqdm import tqdm
import matplotlib.lines as mlines


def get_codes():  # same
    """
    Gets the PheWAS codes from a local csv file and load it into a pandas DataFrame.

    :returns: All of the codes from the resource file.
    :rtype: pandas DataFrame

    """
    sep = os.sep
    path = os.path.dirname(os.path.abspath(__file__))
    filename = os.sep.join([path, 'resources', 'prowas_codes.csv'])
    return pd.read_csv(filename,dtype=str)


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


def get_input(path, filename, reg_type):  # diff -done - add duration
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
    cptfile = pd.read_csv(wholefname)
    cptfile['cpt'] = cptfile['cpt'].str.strip()
    if reg_type == 0:
        phenotypes = pd.merge(cptfile, codes, on='cpt')
        phenotypes['MaxAgeAtCPT'] = 0
        phenotypes['MaxAgeAtCPT'] = phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('max')
    else:
        """
        This needs to be changed, need to adjust for a variety of different naming conventions
        in the phenotype file, not simply 'AgeAtCPT', 'id', 'cpt', etc.
        Either we need to adjust for different names in the code, or state explicitly in the
        documentation that we cannot do things like this.
        """
        phenotypes = pd.merge(cptfile, codes, on='cpt')
        phenotypes['count'] = 0
        phenotypes['count'] = phenotypes.groupby(['id', 'prowas_code'])['count'].transform('count')
        phenotypes['duration'] = phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('max') - \
                                 phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('min') + 1
        phenotypes['MaxAgeAtCPT'] = 0
        phenotypes['MaxAgeAtCPT'] = phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('max')
    return phenotypes


def generate_feature_matrix(genotypes, phenotypes, reg_type, phewas_cov=''):  # diff - done
    """
    Generates the feature matrix that will be used to run the regressions.

    :param genotypes:
    :param phenotypes:
    :type genotypes:
    :type phenotypes:

    :returns:
    :rtype:

    """
    pu=phenotypes[['id','prowas_code']].drop_duplicates()
    temp = pd.DataFrame(np.log2(pu['id'].drop_duplicates().count()/pu.groupby('prowas_code')['id'].count()).reset_index())
    temp.rename(columns={'id': 'idf'}, inplace=True)
    prowas_codes2 = pd.merge(prowas_codes, temp, on='prowas_code', how='left')

    feature_matrix = np.zeros((3, genotypes.shape[0], prowas_codes.shape[0]), dtype=int)
    count = 0;
    for i in tqdm(genotypes['id']):
        if reg_type == 0:
            temp = pd.DataFrame(phenotypes[phenotypes['id'] == i][['prowas_code', 'MaxAgeAtCPT']]).drop_duplicates()
            match = prowas_codes2['prowas_code'].isin(list(phenotypes[phenotypes['id'] == i]['prowas_code']))
            feature_matrix[0][count, match[match == True].index] = 1
            age = pd.merge(prowas_codes2, temp, on='prowas_code', how='left')['MaxAgeAtCPT']
            age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeAtVisit']
            assert np.all(np.isfinite(age)), "make sure MaxAgeAtVisit is filled"
            feature_matrix[1][count, :] = age
            if phewas_cov:
                feature_matrix[2][count, :] = int(phewas_cov in list(phenotypes[phenotypes['id'] == i]['prowas_code']))

        else:
            if reg_type == 1:
                temp = pd.DataFrame(
                    phenotypes[phenotypes['id'] == i][['prowas_code', 'MaxAgeAtCPT', 'count']]).drop_duplicates()
                cts = pd.merge(prowas_codes, temp, on='prowas_code', how='left')['count']
                cts[np.isnan(cts)] = 0
                feature_matrix[0][count, :] = cts
                age = pd.merge(prowas_codes2, temp, on='prowas_code', how='left')['MaxAgeAtCPT']
                age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeAtVisit']
                assert np.all(np.isfinite(age)), "make sure MaxAgeAtVisit is filled"
                feature_matrix[1][count, :] = age
                if phewas_cov:
                    feature_matrix[2][count, :] = int(
                        phewas_cov in list(phenotypes[phenotypes['id'] == i]['prowas_code']))

            elif reg_type == 2:
                temp = pd.DataFrame(
                    phenotypes[phenotypes['id'] == i][['prowas_code', 'MaxAgeAtCPT', 'count']]).drop_duplicates()
                temp = pd.merge(prowas_codes2, temp, on='prowas_code', how='left')
                tfidf=temp['count']*temp['idf']
                tfidf[np.isnan(tfidf)] = 0

                feature_matrix[0][count, :] = tfidf
                age = pd.merge(prowas_codes2, temp, on='prowas_code', how='left')['MaxAgeAtCPT']
                age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeAtVisit']
                assert np.all(np.isfinite(age)), "make sure MaxAgeAtVisit is filled"
                feature_matrix[1][count, :] = age
                if phewas_cov:
                    feature_matrix[2][count, :] = int(
                        phewas_cov in list(phenotypes[phenotypes['id'] == i]['prowas_code']))

        count += 1
    return feature_matrix


"""

Statistical Modeling

"""


def get_phewas_info(p_index):  # same
    """
    Returns all of the info of the phewas code at the given index.

    :param p_index: The index of the desired phewas code
    :type p_index: int

    :returns: A list including the code, the name, and the rollup of the phewas code. The rollup is a list of all of the cpt-9 codes that are grouped into this phewas code.
    :rtype: list of strings
    """
    p_code = prowas_codes.loc[p_index].prowas_code
    corresponding = codes[codes.prowas_code == p_code]

    p_name = corresponding.iloc[0].prowas_desc
    p_rollup = ','.join(codes[codes.prowas_code == p_code].cpt.tolist())
    return [p_code, p_name, p_rollup]


def calculate_odds_ratio(genotypes, phen_vector1, phen_vector2, covariates, lr=0, response='',
						 phen_vector3=''):  # diff - done
	"""
	Runs the regression for a specific phenotype vector relative to the genotype data and covariates.

	:param genotypes: a DataFrame containing the genotype information
	:param phen_vector1: a array containing the phenotype vector
	:param phen_vector2: a array containing the phenotype vector
	:param covariates: a string containing all desired covariates
	:type genotypes: pandas DataFrame
	:type phen_vector1: numpy array
	:type phen_vector2: numpy array
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
	data['MaxAgeAtCPT'] = phen_vector2
	# f='y~'+covariates
	if covariates is not '':
		covariates = '+' + covariates
	if response:
		f = response + '~ y + genotype' + covariates
		if phen_vector3.any():
			data['phe'] = phen_vector3
			f = response + '~ y + phe + genotype' + covariates
	else:
		f = 'genotype ~ y' + covariates
		if phen_vector3.any():
			data['phe'] = phen_vector3
			f = 'genotype ~ y + phe' + covariates
	try:
		if lr == 0: # fit logit without regulatization
			logreg = smf.logit(f, data).fit(disp=False)
			p = logreg.pvalues.y
			odds = 0  #
			conf = logreg.conf_int()
			od = [-math.log10(p), p, logreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
		elif lr == 1: # fit logit with regularization
			f1 = f.split(' ~ ')
			f1[1] = f1[1].replace(" ", "")
			logit = sm.Logit(data[f1[0].strip()], data[f1[1].split('+')])
			lf = logit.fit_regularized(method='l1', alpha=0.1, disp=0, trim_mode='size', qc_verbose=0)
			p = lf.pvalues.y
			odds = 0
			conf = lf.conf_int()
			od = [-math.log10(p), p, lf.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
		else:
			linreg = smf.logit(f, data).fit(method='bfgs', disp=False)
			p = linreg.pvalues.y
			odds = 0
			conf = linreg.conf_int()
			od = [-math.log10(p), p, linreg.params.y, '[%s,%s]' % (conf[0]['y'], conf[1]['y'])]
	except ValueError as ve:
		print(ve)
		print('lr = % d' %lr)
		odds = 0
		p = np.nan
		od = [np.nan, p, np.nan, np.nan]
	except Exception as e:
		print(e)
		odds = 0
		p = np.nan
		od = [np.nan, p, np.nan, np.nan]
	return (odds, p, od)


def run_phewas(fm, genotypes, covariates, reg_type, response='', phewas_cov=''):  # same
	"""
	For each phewas code in the feature matrix, run the specified type of regression and save all of the resulting p-values.

	:param fm: The phewas feature matrix.
	:param genotypes: A pandas DataFrame of the genotype file.
	:param covariates: The covariates that the function is to be run on.
	:param reg_type: The covariates that the function is to be run on.
	:param response: The covariates that the function is to be run on.
	:param phewas_cov: The covariates that the function is to be run on.

	:returns: A tuple containing indices, p-values, and all the regression data.
	"""

	num_phecodes = len(fm[0, 0])
	thresh = math.ceil(genotypes.shape[0] * 0.03)
	# store all of the pertinent data from the regressions
	regressions = pd.DataFrame(columns=output_columns)
	control = fm[0][genotypes.genotype == 0, :]
	disease = fm[0][genotypes.genotype == 1, :]
	# find all phecodes that only present for a single genotype (ie only controls or only diseased show the phecode) -> have to use regularization
	inds = np.where((control.any(axis=0) & ~disease.any(axis=0)) | (~control.any(axis=0) & disease.any(axis=0)))[0]
	for index in tqdm(range(num_phecodes), desc='Running Regressions'):
		phen_vector1 = fm[0][:, index]
		phen_vector2 = fm[1][:, index]
		phen_vector3 = fm[2][:, index]
		# to prevent false positives, only run regressions if more than thresh records have positive values
		if np.where(phen_vector1 > 0)[0].shape[0] > thresh:
			if index in inds:
				res = calculate_odds_ratio(genotypes, phen_vector1, phen_vector2, covariates,
										   lr=1,
										   response=response,
										   phen_vector3=phen_vector3)
			else:
				res = calculate_odds_ratio(genotypes, phen_vector1, phen_vector2, covariates,
										   lr=0,
										   response=response,
										   phen_vector3=phen_vector3)

		else: # default (non-significant) values if not enough samples to run regression
			odds = 0
			p = np.nan
			od = [np.nan, p, np.nan, np.nan]
			res = (odds, p, od)

		# save all of the regression data
		phewas_info = get_phewas_info(index)
		stat_info = res[2]
		info = phewas_info[0:2] + stat_info + [phewas_info[2]]

		regressions.loc[index] = info

	return regressions.dropna(subset=['p-val']).sort_values(by='PheWAS Code')


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
        thresh = 0.05 * i / len(sn)
        if sn[i] <= power:
            break
    return sn[i]


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


def plot_manhattan(regressions, thresh, show_imbalance=True, save='', save_format=''):  # same
	"""
	Plots the data on a Manhattan Plot.

	:param regressions: dataframe containing the regression results
	:param thresh: the significance threshold
	:param save: the output file to save to (if empty, display the plot)
	:param show_imbalance: boolean variable that determines whether or not to show imbalances on the plot (default True)
	:type regressions: pandas DataFrame
	:type thresh: float
	:type save: str
	:type show_imbalance: boolean

	"""

	# Initialize figure
	fig = plt.figure(1)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# Merge regressions with Phewas data to get categories
	regressions = pd.merge(regressions, prowas_codes, left_on='PheWAS Code', right_on='prowas_code').sort_values(by='ccs')

	# Determine whether or not to show the imbalance.
	# show_imbalance = imbalances.size != 0

	# c = icd9_codes.loc[phewas_codes['index']]
	# c = c.reset_index()
	# idx = c.sort_values(by='category').index

	# Plot all points w/ labels
	e = 1
	artists = []
	plt.ylabel('-log10(p)')
	ax.axhline(y=thresh, color='red', ls='dotted')  # plot threshold
	for ix,data in regressions.iterrows():
		logp_ix = data['"-log(p)"']
		if  logp_ix > thresh:
			# determine marker type based on whether/not showing imbalance
			if show_imbalance:
				mew = 1.5
				if data['beta'] > 0: m = '+'
				else: m = '_'
			else:
				mew = 0.0
				m = 'o'
			# Plot PheCode data point & format PheCode label
			ax.plot(e, logp_ix, m, color="blue", fillstyle='full', markeredgewidth=mew)
			artists.append(ax.text(e, logp_ix, data['prowas_desc'], rotation=89, va='bottom', fontsize=6))
			e += 15

	# # Legend
	# line1 = []
	# box = ax.get_position()
	# ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
	# for lab in plot_colors.keys():
	# 	line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=plot_colors[lab], label=lab))
	# artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, 0), loc='upper center', fancybox=True, ncol=4, prop={'size': 6}))

	# Plot x axis
	ax.axhline(y=0, color='black')
	frame1.axes.get_xaxis().set_visible(False)

	# If the imbalance is to be shown, draw lines to show the categories.
	# if show_imbalance:
	# 	for pos in linepos:
	# 		ax.axvline(x=pos, color='black', ls='dotted')

	# Save the plot
	if save:
		plt.savefig(save,format=save_format, bbox_extra_artists=artists, bbox_inches='tight')
		plt.clf()

	return


def plot_odds_ratio(regressions, thresh, show_imbalance=True, save='', save_format='', label_loc="plot"):  # same
	"""
	Plots the data on a Log Odds Plot.

	:param regressions: dataframe containing the regression results
	:param thresh: the significance threshold
	:param save: the output file to save to (if empty, display the plot)
	:param show_imbalance: boolean variable that determines whether or not to show imbalances on the plot (default True)
	:param label_loc: the output file to save to (if empty, display the plot)
	:type regressions: pandas DataFrame
	:type thresh: float
	:type save: str
	:type show_imbalance: boolean

	"""

	# Initialize figure
	fig = plt.figure(2)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# Merge regressions with Phewas data to get categories
	regressions = pd.merge(regressions, prowas_codes, left_on='PheWAS Code', right_on='prowas_code').sort_values(by='ccs')

	# determine whether or not to show imbalances
	# show_imbalance = imbalances.size != 0

	# Sort the phewas codes by category.
	# c = icd9_codes.loc[phewas_codes['index']]
	# c = c.reset_index()
	# idx = c.sort_values(by='category').index

	# Plot all points w/ labels
	e = 1  # vertical index
	ho = 0.025  # horizontal text offset
	vo = 1  # vertical text offset
	text_size = 6
	artists = []
	if label_loc == "axis":
		phecode_labels = []
		phecode_locs = []
	plt.xlabel('Log odds ratio')
	for ix, data in regressions.iterrows():
		beta_ix = data['beta']
		if data['"-log(p)"'] > thresh:
			# Add Phecode label
			if label_loc == "plot":
				if show_imbalance:
					if beta_ix > 0:
						artists.append(
							ax.text(beta_ix + ho, e + vo, data['prowas_desc'], rotation=0, ha='left', fontsize=text_size))
					else:
						artists.append(ax.text(beta_ix - ho, e + vo, data['prowas_desc'], rotation=0, ha='right',
						                       fontsize=text_size))
				else:
					artists.append(
						ax.text(beta_ix + ho, e + vo, data['prowas_desc'], rotation=0, va='bottom', fontsize=text_size))
			else:  # location = "axis"
				phecode_labels.append(data['prowas_desc'])
				phecode_locs.append(e)

			# Plot Phecode Data
			ax.plot(beta_ix, e, 'o', color="green", fillstyle='full', markeredgewidth=0.0)
			ax.plot([data['lowlim'], data['uplim']], [e, e], color="green")
			e += 15

	# Plot y axis
	ax.axvline(x=0, color='black')

	if label_loc == "axis":
		plt.yticks(phecode_locs, phecode_labels, ha='right', fontsize=text_size)
	else:
		frame1.axes.get_yaxis().set_visible(False)

	# Legend
	# line1 = []
	# box = ax.get_position()
	# ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
	# for lab in plot_colors.keys():
	# 	line1.append(
	# 		mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=plot_colors[lab], label=lab))
	# artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, -0.125), loc='upper center', fancybox=True, ncol=4,
	#                          prop={'size': text_size}))

	# If the imbalance is to be shown, draw lines to show the categories.
	# if show_imbalance:
	# 	for pos in linepos:
	# 		ax.axvline(x=pos, color='black', ls='dotted')

	# Save the plot
	if save:
		plt.savefig(save, format=save_format, bbox_extra_artists=artists, bbox_inches='tight')
		plt.clf()

	return

def process_args(kwargs, optargs, *args):
    clean = np.vectorize(lambda x: x[x.rfind('-') + 1:] + '=')
    searchfor = clean(list(optargs.keys()))
    opts, rem = getopt.getopt(args, '', searchfor)
    assert len(rem) == 0, 'Unknown arguments included %s' % (str(rem))
    for option in opts:
        k, v = option
        kwargs[optargs[k]] = v

    return kwargs


def display_kwargs(kwargs):
    print ("Arguments: ")
    for k, v in kwargs.items():
        left = str(k).ljust(30, '.')
        right = str(v).rjust(50, '.')
        print(left + right)


output_columns = ['PheWAS Code',
                  'PheWAS Name',
                  'p-val',
                  '\"-log(p)\"',
                  'beta',
                  'Conf-interval beta',
                  'cpt']

plot_colors = {'-': 'gold',
               'circulatory system': 'red',
               'congenital anomalies': 'mediumspringgreen',
               'dermatologic': 'maroon',
               'digestive': 'green',
               'endocrine/metabolic': 'darkred',
               'genitourinary': 'black',
               'hematopoietic': 'orange',
               'infectious diseases': 'blue',
               'injuries & poisonings': 'slategray',
               'mental disorders': 'fuchsia',
               'musculoskeletal': 'darkgreen',
               'neoplasms': 'teal',
               'neurological': 'midnightblue',
               'pregnancy complications': 'gold',
               'respiratory': 'brown',
               'sense organs': 'darkviolet',
               'symptoms': 'darkviolet'}
imbalance_colors = {
    0: 'white',
    1: 'deepskyblue',
    -1: 'red'
}
regression_map = {
    'log': 0,
    'lin': 1,
    'lind': 2
}
threshold_map = {
    'bon': 0,
    'fdr': 1
}
global codes, prowas_codes
codes = get_codes()
prowas_codes = codes[['prowas_code','prowas_desc','ccs','CCS Label']].drop_duplicates(subset='prowas_code')
prowas_codes.reset_index(inplace=True,drop=True)


