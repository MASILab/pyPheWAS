"""
PyPheWAS Core version 2 (main PyPheWAS code)
Developed by:
    Shikha Chaganti, PhD
    Cailey Kerley

MASI Lab
Department of Electrical Engineering and Computer Science
Vanderbilt University
"""

from collections import Counter
import getopt
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import pandas as pd
import statsmodels.discrete.discrete_model as sm
import statsmodels.formula.api as smf
import matplotlib.lines as mlines
from tqdm import tqdm
import sys


def get_codes(version):  # same
	"""
	Gets the PheWAS codes from a local csv file and load it into a pandas DataFrame.

	:returns: All of the codes from the resource file.
	:rtype: pandas DataFrame

	"""
	path = os.path.dirname(os.path.abspath(__file__))
	filename = os.sep.join([path, 'resources', 'phecode_map_'+ version + '.csv'])
	if '9' in version:
		icd_col = 'ICD9'
	elif '10' in version:
		icd_col = 'ICD10'
	else:
		# TODO: add an actual exception here
		print('Error in phecode map version name: exiting pyPheWAS')
		sys.exit()
	try:
		phecode_map = pd.read_csv(filename, dtype={icd_col: str, 'PheCode': str})
		phecode_map.rename(columns={icd_col:'ICD_CODE'},inplace=True)
		phecode_map.dropna(subset=['PheCode'], inplace=True)
	except Exception as e:
		print(e)
		print('Error in phecode map version name: exiting pyPheWAS')
		sys.exit()

	return phecode_map


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
	genotypes_df = pd.read_csv(wholefname)
	genotypes_df = genotypes_df.dropna(subset=['id'])
	genotypes_df.sort_values(by='id', inplace=True)
	return genotypes_df


def get_icd_codes(path, filename, reg_type):
	"""
	Read all of the phenotype data from the given file and load it into a pandas DataFrame.

	:param path: The path to the file that contains the phenotype data
	:param filename: The name of the file that contains the phenotype data.
	:param reg_type: Type of regression
	:type path: string
	:type filename: string
	:type reg_type: int

	:returns: The data from the phenotype file.
	:rtype: pandas DataFrame
	"""

	wholefname = path + filename
	icdfile = pd.read_csv(wholefname,dtype={'ICD_CODE':str})
	icdfile['ICD_CODE'] = icdfile['ICD_CODE'].str.strip()
	icd_types = np.unique(icdfile['ICD_TYPE'])

	# check ICD types present in file
	if not all((icd_types == 9) | (icd_types == 10)):
		# TODO: add actual exception
		print('Found an ICD_TYPE that was not 9 or 10 - Please check phenotype file.\nExiting pyPheWAS')
		sys.exit()

	# merge with Phecode table, depends on which type(s) of ICD codes are in the icd file
	if icd_types.shape[0] == 2:
		print('Found both ICD-9 and ICD-10 codes.')
		icd9s = icdfile[icdfile['ICD_TYPE'] == 9]
		icd10s = icdfile[icdfile['ICD_TYPE'] == 10]
		phenotypes_9 = pd.merge(icd9s, icd9_codes,on='ICD_CODE')
		phenotypes_10 = pd.merge(icd10s, icd10_codes, on='ICD_CODE')
		phenotypes = phenotypes_9.append(phenotypes_10,sort=False)
	elif icd_types[0] == 9:
		print('Found only ICD-9 codes.')
		phenotypes = pd.merge(icdfile,icd9_codes,on='ICD_CODE')
	elif icd_types[0] == 10:
		print('Found only ICD-10 codes.')
		phenotypes = pd.merge(icdfile, icd10_codes, on='ICD_CODE')
	else:
		# TODO: add actual exception
		print('An issue occurred while parsing the ICD_TYPE column - Please check phenotype file.\nExiting pyPheWAS')
		sys.exit()

	if reg_type == 2:
		# pre-calculate durations for feature matrix step
		phenotypes['duration'] = phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('max') - \
								 phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('min') + 1

	# calculate max age at each ICD code
	phenotypes['MaxAgeAtICD'] = 0
	phenotypes['MaxAgeAtICD'] = phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('max')
	# pre-sort by ID for feature matrix step
	phenotypes.sort_values(by='id', inplace=True)

	return phenotypes


def generate_feature_matrix(genotypes_df, icds, reg_type, phewas_cov):
	"""
	Generates the feature matrix that will be used to run the regressions.

	:param genotypes_df:
	:param icds:
	:param reg_type:
	:param phewas_cov:
	:type genotypes_df: pandas DataFrame
	:type icds: pandas DataFrame
	:type reg_type: int
	:type phewas_cov: str

	:returns:
	:rtype:

	"""

	# feature matrix 0: reg_type specific data
	# feature matrix 1: age data
	# feature matrix 2: phecode covariance data
	feature_matrix = np.zeros((3, genotypes_df.shape[0], phewas_codes.shape[0]), dtype=float)

	# make genotype a dictionary for faster access time
	genotypes = genotypes_df.set_index('id').to_dict('index')

	# use phewascodes to make a dictionary of indices in the np array
	empty_phewas_df = phewas_codes.set_index('PheCode')
	empty_phewas_df.sort_index(inplace=True)
	empty_phewas_df['np_index'] = range(0,empty_phewas_df.shape[0])
	np_index = empty_phewas_df['np_index'].to_dict()

	exclude = []  # list of ids to exclude (in icd list but not in genotype list)
	last_id = ''  # track last id seen in icd list
	count = -1

	for _,event in tqdm(icds.iterrows(), desc="Processing ICDs", total=icds.shape[0]):
		curr_id = event['id']
		if not curr_id in genotypes:
			if not curr_id in exclude:
				print('%s has records in icd file but is not in group file - excluding from study' % curr_id)
				exclude.append(curr_id)
			continue
		# check id to see if a new subject has been found
		if last_id != curr_id:
			count += 1
			last_id = curr_id  # reset last_id
			feature_matrix[1][count] = genotypes[curr_id]['MaxAgeAtVisit']

		# get column index of event phecode
		phecode_ix = np_index[event['PheCode']]
		# add the max at of that ICD code to the age feature matrix
		feature_matrix[1][count][phecode_ix] = event['MaxAgeAtICD']
		# if using a phecode covariate and the event phecode is equal to the covariate phecode, set the fm2 to 1 for this subject
		if (phewas_cov is not None) & (event['PheCode'] == phewas_cov):
			feature_matrix[2][count][:] = 1

		if reg_type == 0:
			# binary aggregate: add a 1 to the phecode's column to denote that this subject had this phecode at some point
			feature_matrix[0][count][phecode_ix] = 1
		elif reg_type == 1:
			# count aggregate: add 1 to the phecode's column to find the total number of times this subject had this phecode
			feature_matrix[0][count][phecode_ix] += 1
		else:
			# dur aggregate: store the number of years between the first and last events with this phecode
			feature_matrix[0][count][phecode_ix] = event['duration'] # duration calculated in get_icd_codes()

	return feature_matrix



"""

Statistical Modeling

"""


def get_phewas_info(p_index):  # same
	"""
	Returns all of the info of the phewas code at the given index.

	:param p_index: The index of the desired phewas code
	:type p_index: int

	:returns: A list including the code, the name, and the rollup of the phewas code. The rollup is a list of all of the ICD-9 and ICD-10 codes that are grouped into this phewas code.
	:rtype: list of strings
	"""
	p_code = phewas_codes.loc[p_index, 'PheCode']
	p_name = phewas_codes.loc[p_index, 'Phenotype']
	icd9_ix = icd9_codes['PheCode'] == p_code
	icd10_ix = icd10_codes['PheCode'] == p_code

	p_icd9_rollup = '/'.join(icd9_codes.loc[icd9_ix, 'ICD_CODE'].tolist())
	p_icd10_rollup = '/'.join(icd10_codes.loc[icd10_ix, 'ICD_CODE'].tolist())

	return [p_code, p_name, p_icd9_rollup,p_icd10_rollup]


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
	data['MaxAgeAtICD'] = phen_vector2
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
		od = [np.nan, np.nan, np.nan, np.nan]
	except Exception as e:	
		print(e)
		odds = 0
		p = np.nan
		od = [np.nan, np.nan, np.nan, np.nan]
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
	thresh = math.ceil(genotypes.shape[0] * 0.01)
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
			p = 1
			od = [-0.0, 1.0, 0.0, np.nan]
			res = (odds, p, od)

		# save all of the regression data

		phewas_info = get_phewas_info(index)
		stat_info = res[2]
		info = phewas_info[0:2] + stat_info + [phewas_info[2]] + [phewas_info[3]]

		regressions.loc[index] = info

	return regressions


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



def get_bhy_thresh(p_values, power): # TODO: why does this exist
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
		thresh = power * i / (8.1 * len(sn))
		if sn[i] <= thresh:
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


def plot_data_points(y, thresh, save='', imbalances=np.array([])):  # same
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

	# Initialize figure
	fig = plt.figure(1)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# Determine whether or not to show the imbalance.
	show_imbalance = imbalances.size != 0

	# Sort the phewas codes by category.
	c = icd9_codes.loc[phewas_codes['index']]
	c = c.reset_index()
	idx = c.sort_values(by='category').index

	# Plot all points w/ labels
	e = 1
	artists = []
	plt.ylabel('-log10(p)')
	ax.axhline(y=thresh, color='red', ls='dotted')  # plot threshold
	for i in idx:
		if y[i] > thresh:
			# determine marker type based on whether/not showing imbalance
			if show_imbalance:
				mew = 1.5
				if imbalances[i] > 0: m = '+'
				else: m = '_'
			else:
				mew = 0.0
				m = 'o'
			# Plot PheCode data point & format PheCode label
			ax.plot(e, y[i], m, color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full', markeredgewidth=mew)
			artists.append(ax.text(e, y[i], c['Phenotype'][i], rotation=89, va='bottom', fontsize=6))
			e += 15

	# Legend
	line1 = []
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
	for lab in plot_colors.keys():
		line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=plot_colors[lab], label=lab))
	artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, 0), loc='upper center', fancybox=True, ncol=4, prop={'size': 6}))

	# Plot x axis
	ax.axhline(y=0, color='black')
	frame1.axes.get_xaxis().set_visible(False)

	# If the imbalance is to be shown, draw lines to show the categories.
	# if show_imbalance:
	# 	for pos in linepos:
	# 		ax.axvline(x=pos, color='black', ls='dotted')

	# Save the plot
	if save:
		pdf = PdfPages(save)
		pdf.savefig(bbox_extra_artists=artists, bbox_inches='tight')
		pdf.close()
		plt.clf()

	return


def plot_odds_ratio(y, p, thresh, save='', imbalances=np.array([])):  # same
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

	# Initialize figure
	fig = plt.figure(2)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# determine whether or not to show imbalances
	show_imbalance = imbalances.size != 0

	# Sort the phewas codes by category.
	c = icd9_codes.loc[phewas_codes['index']]
	c = c.reset_index()
	idx = c.sort_values(by='category').index

	# Plot all points w/ labels
	e = 1
	artists = []
	plt.xlabel('Log odds ratio')
	for i in idx:
		if p[i] > thresh:
			# Add Phecode label
			if show_imbalance:
				if imbalances[i] > 0:
					artists.append(ax.text(y[i][0], e, c['Phenotype'][i], rotation=0, ha='left', fontsize=6))
				else:
					artists.append(ax.text(y[i][0], e, c['Phenotype'][i], rotation=0, ha='right', fontsize=6))
			else:
				artists.append(ax.text(y[i][0], e, c['Phenotype'][i], rotation=0, va='bottom', fontsize=6))

			# Plot Phecode Data
			ax.plot(y[i][0], e, 'o', color=plot_colors[c[i:i + 1].category_string.values[0]], fillstyle='full', markeredgewidth=0.0)
			ax.plot([y[i, 1], y[i, 2]], [e, e], color=plot_colors[c[i:i + 1].category_string.values[0]])
			e += 15

	# Legend
	line1 = []
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
	for lab in plot_colors.keys():
		line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=plot_colors[lab], label=lab))
	artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, -0.125), loc='upper center', fancybox=True, ncol=4, prop={'size': 6}))

	# Plot y axis
	ax.axvline(x=0, color='black')
	frame1.axes.get_yaxis().set_visible(False)

	# If the imbalance is to be shown, draw lines to show the categories.
	# if show_imbalance:
	# 	for pos in linepos:
	# 		ax.axvline(x=pos, color='black', ls='dotted')

	# Save the plot
	if save:
		pdf = PdfPages(save)
		pdf.savefig(bbox_extra_artists=artists, bbox_inches='tight')
		pdf.close()
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
	print("Arguments: ")
	for k, v in kwargs.items():
		num_dots = 80 - len(str(k)) - len(str(v))
		# left = str(k).ljust(30, '.')
		# right = str(v).rjust(50, '.')
		print(str(k) + '.'*num_dots + str(v))


output_columns = ['PheWAS Code',
				  'PheWAS Name',
				  '\"-log(p)\"',
				  'p-val',
				  'beta',
				  'Conf-interval beta',
				  'ICD-9',
				  'ICD-10']

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
	'dur': 2
}
threshold_map = {
	'bon': 0,
	'fdr': 1,
	'custom':2
}

icd9_codes = get_codes(version='v1_2_icd9')
icd10_codes = get_codes(version='v1_2_icd10_beta')
phewas_codes = icd9_codes[['PheCode','Phenotype','category','category_string']].drop_duplicates(subset='PheCode')
phewas_codes.sort_values(by=['PheCode'], inplace=True)
phewas_codes.reset_index(inplace=True)
