"""
PyPheWAS Explorer Core Functions
Developed by:
    Cailey Kerley

MASI Lab
Department of Electrical Engineering and Computer Science
Vanderbilt University
"""

from pyPheWAS.pyPhewasCorev2 import icd9_codes, icd10_codes, phewas_codes
from collections import Counter
import getopt
import math
import numpy as np
import os
import pandas as pd
import statsmodels.discrete.discrete_model as sm
import statsmodels.formula.api as smf
from tqdm import tqdm
import sys
import scipy.stats
import http.server
import socketserver


output_columns = ['PheCode',
				  'Phenotype',
				  'Pheno_id',
				  'count',
				  'neg_log_p',
				  'pval',
				  'beta',
				  'beta_ci_low',
				  'beta_ci_up',
				  'category_id',
				  'category',
				  'ICD9',
				  'ICD10']


"""

Preprocessing

"""

def get_group_file(filename):  # same
	"""
	Read all of the genotype data from the given file and load it into a pandas DataFrame.

	:param path: The path to the file that contains the phenotype data
	:param filename: The name of the file that contains the phenotype data.
	:type path: string
	:type filename: string

	:returns: The data from the genotype file.
	:rtype: pandas DataFrame
	"""
	genotypes_df = pd.read_csv(filename)
	genotypes_df = genotypes_df.dropna(subset=['id'])
	genotypes_df.sort_values(by='id', inplace=True)
	genotypes_df.reset_index(inplace=True, drop=True)
	return genotypes_df


def get_icd_codes(filename):
	"""
	Read all of the phenotype data from the given file and load it into a pandas DataFrame.

	:param filename: The name of the file that contains the phenotype data.
	:type filename: string

	:returns: The data from the phenotype file.
	:rtype: pandas DataFrame
	"""

	icdfile = pd.read_csv(filename,dtype={'ICD_CODE':str})
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

	# pre-calculate durations for feature matrix step
	phenotypes['duration'] = phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('max') - \
							 phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('min') + 1

	phenotypes.sort_values(by='id', inplace=True)

	return phenotypes


def generate_feature_matrix(genotypes, icds):
	"""
	Generates the feature matrix that will be used to run the regressions.

	:param genotypes_df:
	:param icds:
	:type genotypes_df: pandas DataFrame
	:type icds: pandas DataFrame

	:returns:
	:rtype:

	"""
	# init
	fm_bin = np.zeros((genotypes.shape[0], phewas_codes.shape[0]), dtype=float)
	fm_cnt = np.zeros((genotypes.shape[0], phewas_codes.shape[0]), dtype=float)
	fm_dur = np.zeros((genotypes.shape[0], phewas_codes.shape[0]), dtype=float)

	# make genotype a dictionary for faster access time
	genotypes_dict = genotypes.set_index('id').to_dict('index')

	# use phewas codes to make a dictionary of indices in the np array
	empty_phewas_df = phewas_codes.set_index('PheCode')
	empty_phewas_df.sort_index(inplace=True)
	empty_phewas_df['np_index'] = range(0,empty_phewas_df.shape[0])
	np_index = empty_phewas_df['np_index'].to_dict()

	exclude = []  # list of ids to exclude (in icd list but not in genotype list)
	last_id = ''  # track last id seen in icd list
	count = -1

	for _,event in tqdm(icds.iterrows(), desc="Processing ICDs", total=icds.shape[0]):
		curr_id = event['id']
		if not curr_id in genotypes_dict:
			if not curr_id in exclude:
				print('%s has records in icd file but is not in group file - excluding from study' % curr_id)
				exclude.append(curr_id)
			continue
		# check id to see if a new subject has been found
		if last_id != curr_id:
			count += 1
			while(curr_id != genotypes.loc[count,'id']): # subject at genotypes.loc[count] does not have ICD codes
				count += 1 # continue processing next subject
			last_id = curr_id  # reset last_id

		# get column index of event phecode
		phecode_ix = np_index[event['PheCode']]

		# binary aggregate: add a 1 to the phecode's column to denote that this subject had this phecode at some point
		fm_bin[count][phecode_ix] = 1
		# count aggregate: add 1 to the phecode's column to find the total number of times this subject had this phecode
		fm_cnt[count][phecode_ix] += 1
		# dur aggregate: store the number of years between the first and last events with this phecode
		fm_dur[count][phecode_ix] = event['duration']

	return fm_bin, fm_cnt, fm_dur, list(empty_phewas_df.index)


"""
Group Var Analysis
"""

def get_1D_histogram(group, var_name, response):

	bin_max = max(group[var_name])
	bin_min = min(group[var_name])

	geno0 = group[group[response] == 0]
	geno1 = group[group[response] == 1]

	H0, edges0 = np.histogram(geno0[var_name].to_numpy(),range=(bin_min,bin_max))
	H1, edges1 = np.histogram(geno1[var_name].to_numpy(),range=(bin_min,bin_max))

	H_df0 = pd.DataFrame(index=range(0,len(H0)))
	H_df0['count'] = H0
	H_df0['xmin'] = edges0[0:-1]
	H_df0['xmax'] = edges0[1:]
	H_df0['response'] = 0

	H_df1 = pd.DataFrame(index=range(0, len(H1)))
	H_df1['count'] = H1
	H_df1['xmin'] = edges1[0:-1]
	H_df1['xmax'] = edges1[1:]
	H_df1['response'] = 1

	H_df = H_df0.append(H_df1)

	H_df['var_name'] = var_name
	return H_df.reset_index(drop=True)


def get_2D_histogram(group, var1, var2, response):

	v1_bin_max = max(group[var1])
	v1_bin_min = min(group[var1])
	v2_bin_max = max(group[var2])
	v2_bin_min = min(group[var2])

	bin_range = [[v1_bin_min,v1_bin_max],[v2_bin_min,v2_bin_max]]

	geno0 = group[group[response] == 0]
	geno1 = group[group[response] == 1]

	H0, v1_edges0, v2_edges0 = np.histogram2d(geno0[var1].to_numpy(),geno0[var2].to_numpy(),range=bin_range)
	H1, v1_edges1, v2_edges1 = np.histogram2d(geno1[var1].to_numpy(),geno1[var2].to_numpy(),range=bin_range)

	H_df = pd.DataFrame(columns=['response','v1_ix','v2_ix','count','v1_min','v1_max','v2_min','v2_max'])
	k = 0
	for i in range(0, H0.shape[0]):
		for j in range(0, H0.shape[1]):
			# response = 0
			H_df.loc[k, 'response'] = 0
			H_df.loc[k, 'v1_ix'] = i
			H_df.loc[k, 'v2_ix'] = j
			H_df.loc[k, 'count'] = H0[i,j]
			H_df.loc[k, 'v1_min'] = v1_edges0[i]
			H_df.loc[k, 'v1_max'] = v1_edges0[i+1]
			H_df.loc[k, 'v2_min'] = v2_edges0[j]
			H_df.loc[k, 'v2_max'] = v2_edges0[j+1]
			k+=1
			# response = 1
			H_df.loc[k, 'response'] = 1
			H_df.loc[k, 'v1_ix'] = i
			H_df.loc[k, 'v2_ix'] = j
			H_df.loc[k, 'count'] = H1[i,j]
			H_df.loc[k, 'v1_min'] = v1_edges1[i]
			H_df.loc[k, 'v1_max'] = v1_edges1[i+1]
			H_df.loc[k, 'v2_min'] = v2_edges1[j]
			H_df.loc[k, 'v2_max'] = v2_edges1[j+1]
			k+=1

	H_df['v1'] = var1
	H_df['v2'] = var2

	return H_df

def variable_comparison(group, var1, var2, response):
	res_df = pd.DataFrame(columns=['test_name','result','pval','var'])
	g_df = group[[response, var1, var2]].copy()  # copy to avoid modifying original df
	g_df.rename(columns={var1: "var1", var2: "var2"}, inplace=True)
	ix = 0
	# correlation
	res_df.loc[ix, 'test_name'] = 'correlation'
	res_df.loc[ix, 'var'] = 'joint'
	var_data = g_df[["var1", "var2"]].to_numpy()
	[corr, pval] = scipy.stats.spearmanr(var_data)
	res_df.loc[ix, 'result'] = corr
	res_df.loc[ix, 'pval'] = pval
	ix += 1

	# univariate regressions
	res_df.loc[ix, 'test_name'] = 'univariate regression'
	res_df.loc[ix, 'var'] = var1
	f_v1 = response + ' ~ var1'
	logreg = smf.logit(f_v1, g_df).fit(disp=False)
	res_df.loc[ix, 'pval'] = logreg.pvalues.var1
	res_df.loc[ix, 'result'] = logreg.params.var1
	ix+=1
	res_df.loc[ix, 'test_name'] = 'univariate regression'
	res_df.loc[ix, 'var'] = var2
	f_v2 = response + ' ~ var2'
	logreg = smf.logit(f_v2, g_df).fit(disp=False)
	res_df.loc[ix, 'pval'] = logreg.pvalues.var2
	res_df.loc[ix, 'result'] = logreg.params.var2
	ix += 1

	# multivariate regression
	f_mul = response + ' ~ var1 + var2'
	logreg = smf.logit(f_mul, g_df).fit(disp=False)
	res_df.loc[ix, 'test_name'] = 'multivariate regression'
	res_df.loc[ix, 'var'] = var1
	res_df.loc[ix, 'pval'] = logreg.pvalues.var1
	res_df.loc[ix, 'result'] = logreg.params.var1
	ix += 1
	res_df.loc[ix, 'test_name'] = 'multivariate regression'
	res_df.loc[ix, 'var'] = var2
	res_df.loc[ix, 'pval'] = logreg.pvalues.var2
	res_df.loc[ix, 'result'] = logreg.params.var2

	return res_df



"""

Statistical Modeling

"""


def get_phenotype_info(p_index):  # same
	"""
	Returns all of the info of the phewas code at the given index.

	:param p_index: The index of the desired phewas code
	:type p_index: int

	:returns: A list including the code, the name, and the rollup of the phewas code. The rollup is a list of all of the ICD-9 and ICD-10 codes that are grouped into this phewas code.
	:rtype: list of strings
	"""
	p_code = phewas_codes.loc[p_index, 'PheCode']
	p_name = phewas_codes.loc[p_index, 'Phenotype']
	p_id = p_name.replace('(','')
	p_id = p_id.replace(')', '')
	p_id = p_id.replace(',', '')
	p_id = p_id.replace(' ', '')
	p_id = p_id.replace('\'', '')
	p_id = p_id.replace('\"', '')
	p_id = p_id.replace('\\', '')
	p_id = p_id.replace('/', '')

	cat_id = phewas_codes.loc[p_index, 'category']
	cat = phewas_codes.loc[p_index, 'category_string']
	icd9_ix = icd9_codes['PheCode'] == p_code
	icd10_ix = icd10_codes['PheCode'] == p_code

	p_icd9_rollup = '/'.join(icd9_codes.loc[icd9_ix, 'ICD_CODE'].tolist())
	p_icd10_rollup = '/'.join(icd10_codes.loc[icd10_ix, 'ICD_CODE'].tolist())

	return [p_code, p_name, p_id, cat_id, cat, p_icd9_rollup, p_icd10_rollup]


def fit_pheno_model(genotypes, phen_vector, response, covariates='', phenotype=''):
	"""
	Runs the regression for a specific phenotype vector relative to the genotype data and covariates.

	:param genotypes: a DataFrame containing the genotype information
	:param phen_vector: a array containing the aggregate phenotype vector
	:param response: response variable in the logit model
	:param covariates: *[optional]* covariates to include in the regressions separated by '+' (e.g. 'sex+ageAtDx')
	:param phenotype: *[optional]* phenotype info [code, description] for this regression (used only for error handling)
	:type genotypes: pandas DataFrame
	:type phen_vector: numpy array
	:type response: string
	:type covariates: string
	:type phenotype: list of strings

	:returns: regression statistics
	:rtype: list

	.. note::
		The covariates must be a string that is delimited by '+', not a list.
		If you are using a list of covariates and would like to convert it to the pyPhewas format, use the following::

			l = ['genotype', 'age'] # a list of your covariates
			covariates = '+'.join(l) # pyPhewas format

		The covariates that are listed here *must* be headers to your genotype CSV file.
	"""

	data = genotypes.copy()
	data['y'] = phen_vector

	# append '+' to covariates (if there are any) -> makes definition of 'f' more elegant
	if covariates != '':
		covariates = '+' + covariates

	# define model ('f') for the logisitc regression
	predictors = covariates.replace(" ", "").split('+')
	predictors[0] = 'y'
	f = [response.strip(), predictors]

	try:
		# fit logit with regularization
		logit = sm.Logit(data[f[0]], data[f[1]])
		model = logit.fit_regularized(method='l1', alpha=0.1, disp=0, trim_mode='size', qc_verbose=0)
		# get results
		p = model.pvalues.y
		beta = model.params.y
		conf = model.conf_int()
		reg_result = [-math.log10(p), p, beta, conf[0]['y'], conf[1]['y']]  # collect results
	except Exception as e:
		print('\n')
		if phenotype != '':
			print('ERROR computing regression for phenotype %s (%s)' %(phenotype[0],phenotype[1]))
		print(e)
		reg_result = [np.nan, np.nan, np.nan, np.nan, np.nan]
	return reg_result


def run_phewas(fm, genotypes, covariates, response):
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

	num_phecodes = fm.shape[1]
	thresh = 5 # math.ceil(genotypes.shape[0] * 0.05)
	# store all of the pertinent data from the regressions
	regressions = pd.DataFrame(columns=output_columns)
	control = fm[genotypes[response] == 0, :]
	disease = fm[genotypes[response] == 1, :]

	for index in tqdm(range(num_phecodes), desc='Running Regressions'):
		phen_info = get_phenotype_info(index)
		phen_vector = fm[:, index]
		# to prevent false positives, only run regressions if more than thresh records have positive values
		if np.where(phen_vector > 0)[0].shape[0] > thresh:
			stat_info = fit_pheno_model(genotypes, phen_vector, response, covariates=covariates, phenotype=phen_info[0:2])
		else:
			# not enough samples to run regression
			stat_info = [np.nan, np.nan, np.nan, np.nan, np.nan]
		# save all of the regression data
		info = phen_info[0:3] + [np.where(phen_vector > 0)[0].shape[0]] + stat_info + phen_info[3:5] + [phen_info[5]] + [phen_info[6]]
		regressions.loc[index] = info

	return regressions.dropna(subset=['pval']).sort_values(by='PheCode')


"""
Start the JavaScript GUI
"""

def run_Explorer_GUI():
	PORT = 8000

	web_dir = os.path.join(os.path.dirname(__file__), 'Explorer_GUI')

	os.chdir(web_dir)

	Handler = http.server.SimpleHTTPRequestHandler
	httpd = socketserver.TCPServer(("", PORT), Handler)

	httpd.serve_forever()
	return # kind of unnecessay but whatever
