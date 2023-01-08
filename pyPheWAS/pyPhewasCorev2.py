"""
**pyPheWAS Core version 2 (main pyPheWAS code)**

Contains all functions that drive the core PheWAS & ProWAS analysis tools.
"""

from collections import Counter
import getopt
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import statsmodels.discrete.discrete_model as sm
import statsmodels.formula.api as smf
import matplotlib.lines as mlines
from tqdm import tqdm
import sys
import matplotlib
import warnings


def print_start_msg():
	path = os.path.dirname(os.path.abspath(__file__))
	filename = os.sep.join([path, 'resources', 'pyPheWAS_start_msg.txt'])
	with open(filename, 'r') as msg_f:
		print('\n')
		for line in msg_f:
			print(line.strip('\n'))
	return

def get_codes(filename):
	"""
	Load PheWAS/ProWAS code map from the resources directory.

	:param filename: name of file in the resources folder to load
	:type filename: str

	:returns: code map from the given file
	:rtype: pandas DataFrame

	"""
	path = os.path.dirname(os.path.abspath(__file__))
	filename = os.sep.join([path, 'resources', filename])
	if 'icd9' in filename:
		data_col = 'ICD9'
		new_data_col = 'ICD_CODE'
		phecode_col = 'PheCode'
	elif 'icd10' in filename:
		data_col = 'ICD10'
		new_data_col = 'ICD_CODE'
		phecode_col = 'PheCode'
	else:
		data_col = 'cpt'
		new_data_col = 'CPT_CODE'
		phecode_col = 'prowas_code'
	try:
		phecode_map = pd.read_csv(filename, dtype={data_col:str, phecode_col:str})
		phecode_map.rename(columns={data_col:new_data_col},inplace=True)
		phecode_map.dropna(subset=[phecode_col], inplace=True)
		phecode_map.drop_duplicates(subset=[new_data_col,phecode_col], inplace=True)
	except Exception as e:
		print(e.args[0])
		print('Error loading phecode map : exiting pyPheWAS')
		sys.exit()

	return phecode_map


def get_group_file(path, filename):
	"""
	Read group data from the given file.
	Note: Any records with a null **id** are dropped.

	:param path: path to the file that contains the group data
	:param filename: name of the file that contains the group data.
	:type path: pathlib Path
	:type filename: string

	:returns: The data from the group file.
	:rtype: pandas DataFrame
	"""
	wholefname = path / filename
	genotypes = pd.read_csv(wholefname)
	genotypes = genotypes.dropna(subset=['id'])
	return genotypes


def get_icd_codes(path, filename, reg_type):
	"""
	Read ICD data from the given file and load it into a pandas DataFrame.

	ICD records are mapped to their correpsonding PheWAS Codes.
	The maximum age of each subject at each PheWAS Code is calculated and
	added to the DataFrame as the column *MaxAgeAtICD*. If ``reg_type`` = 2, the
	interval of time (years) over which a subject experiences each PheWAS Code
	is added as the column *duration*.

	:param path: path to the file that contains the phenotype data
	:param filename: name of the file that contains the phenotype data.
	:param reg_type: type of regression (0:binary, 1:count, 2:duration)
	:type path: pathlib Path
	:type filename: str
	:type reg_type: int

	:returns: data from the phenotype file.
	:rtype: pandas DataFrame

	"""

	wholefname = path / filename
	icdfile = pd.read_csv(wholefname,dtype={'ICD_CODE':str})
	icdfile['ICD_CODE'] = icdfile['ICD_CODE'].str.strip()
	icd_types = np.unique(icdfile['ICD_TYPE'])

	# check ICD types present in file
	if not all((icd_types == 9) | (icd_types == 10)):
		raise Exception('Found an ICD_TYPE that was not 9 or 10 - Please check phenotype file.')
	
	# merge with Phecode table, depends on which type(s) of ICD codes are in the icd file
	if icd_types.shape[0] == 2:
		print('Found both ICD-9 and ICD-10 codes.')
		icd9s = icdfile[icdfile['ICD_TYPE'] == 9]
		icd10s = icdfile[icdfile['ICD_TYPE'] == 10]
		phenotypes_9 = pd.merge(icd9s, icd9_codes[['ICD_CODE', 'PheCode']],on='ICD_CODE',how='left')
		phenotypes_10 = pd.merge(icd10s, icd10_codes[['ICD_CODE', 'PheCode']], on='ICD_CODE',how='left')
		phenotypes = pd.concat([phenotypes_9, phenotypes_10], sort=False, ignore_index=True)
	elif icd_types[0] == 9:
		print('Found only ICD-9 codes.')
		phenotypes = pd.merge(icdfile,icd9_codes[['ICD_CODE', 'PheCode']],on='ICD_CODE',how='left')
	elif icd_types[0] == 10:
		print('Found only ICD-10 codes.')
		phenotypes = pd.merge(icdfile, icd10_codes[['ICD_CODE', 'PheCode']], on='ICD_CODE',how='left')
	else:
		raise Exception('An issue occurred while parsing the ICD_TYPE column - Please check phenotype file.')

	# check to see if anything was dropped because of incomplete mapping
	na_mask = phenotypes['PheCode'].isna()
	if np.any(na_mask):
		print('WARNING: Some ICD events were dropped during the PheCode mapping process.')
		num_events = icdfile.shape[0]
		dropped = phenotypes[na_mask].copy()
		percent_drop = 100.0 * (float(dropped.shape[0]) / float(num_events))
		print('\t-- Dropped %d out of %d ICD events (%0.3f%%)' % (dropped.shape[0], num_events, percent_drop))
		print('\t-- Saving dropped events to %s' % (path/'dropped_events.csv'))
		dropped.to_csv(path/'dropped_events.csv',index=False)
		phenotypes.dropna(subset=['PheCode'], inplace=True) # actually remove the events

	if reg_type == 2:
		# pre-calculate durations for feature matrix step
		phenotypes['duration'] = phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('max') - \
								 phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('min') + 1

	# calculate max age at each ICD code
	phenotypes['MaxAgeAtICD'] = 0
	phenotypes['MaxAgeAtICD'] = phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('max')

	return phenotypes


def get_cpt_codes(path, filename, reg_type):
	"""
	Read CPT data from the given file and load it into a pandas DataFrame.

	CPT records are mapped to their correpsonding ProWAS Codes.
	The maximum age of each subject at each ProWAS Code is calculated and
	added to the DataFrame as the column *MaxAgeAtCPT*. If ``reg_type`` = 2, the
	interval of time (years) over which a subject experiences each ProWAS Code
	is added as the column *duration*.

	:param path: path to the file that contains the phenotype data
	:param filename: name of the file that contains the phenotype data.
	:param reg_type: type of regression (0:binary, 1:count, 2:duration)
	:type path: pathlib Path
	:type filename: str
	:type reg_type: int

	:returns: data from the phenotype file.
	:rtype: pandas DataFrame

	"""

	wholefname = path / filename
	cptfile = pd.read_csv(wholefname,dtype={'CPT_CODE':str})
	cptfile['CPT_CODE'] = cptfile['CPT_CODE'].str.strip()
	phenotypes = pd.merge(cptfile, cpt_codes[['CPT_CODE','prowas_code']], on='CPT_CODE', how='left')

	# check to see if anything was dropped because of incomplete mapping
	na_mask = phenotypes['prowas_code'].isna()
	if np.any(na_mask):
		print('WARNING: Some CPT events were dropped during the ProCode mapping process.')
		num_events = cptfile.shape[0]
		dropped = phenotypes[na_mask].copy()
		percent_drop = 100.0 * (float(dropped.shape[0]) / float(num_events))
		print('\t-- Dropped %d out of %d CPT events (%0.3f%%)' % (dropped.shape[0], num_events, percent_drop))
		print('\t-- Saving dropped events to %s' % (path/'dropped_events.csv'))
		dropped.to_csv(path/'dropped_events.csv',index=False)
		phenotypes.dropna(subset=['prowas_code'], inplace=True) # actually remove the events

	if reg_type == 2:
		# pre-calculate durations for feature matrix step
		phenotypes['duration'] = phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('max') - \
								 phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('min') + 1

	# calculate max age at each ICD code
	phenotypes['MaxAgeAtCPT'] = 0
	phenotypes['MaxAgeAtCPT'] = phenotypes.groupby(['id', 'prowas_code'])['AgeAtCPT'].transform('max')

	return phenotypes


def generate_feature_matrix(genotypes, phenotype, reg_type, code_type, pheno_cov=None):
	"""
	Generates the feature matrix that will be used to run the PheWAS or ProWAS analysis. Feature matrix is 3xNxP,
	where N = number of subjects and P = number of PheWAS/ProWAS Codes.

	* Feature matrix [0] is the aggregate PheWAS/ProWAS matrix. It contains phenotype data aggregated across
	  each subject's record accoring to the specified ``reg_type``.
	* Feature matrix [1] is the age matrix, containing each subject's maximum age recorded for each phenotype.
	  For phenotypes absent in the subject's record, the subject's overall maximum age is recorded.
	* Feature matrix [2] is the phenotype covariate matrix. If ``pheno_cov`` is defined, it records whether or
	  not each subject has at least one instance of the specified phenotype in their record. Otherwise, it's a zero matrix.


	:param genotypes: group data
	:param phenotype: phenotype data retrieved from ``pyPheWAS.pyPhewasCorev2.get_icd_codes`` or ``pyPheWAS.pyPhewasCorev2.get_cpt_codes``
	:param reg_type: type of regression (0:binary, 1:count, 2:duration)
	:param code_type: type of EMR code ('ICD' or 'CPT')
	:param pheno_cov: *[optional]* PheWAS or ProWAS code to use as a covariate
	:type genotypes: pandas DataFrame
	:type phenotype: pandas DataFrame
	:type reg_type: int
	:type code_type: str
	:type pheno_cov: str

	:returns (feature_matrix, phenotype_header): (3xNxP feature matrix, PheWAS/ProWAS codes that correspond to the columns in the feature matrix)
	:rtype: (numpy array, list of str)

	.. note:: If ICD/CPT codes do not have a corresponding PheWAS/ProWAS code, they are removed by
           ``pyPheWAS.pyPhewasCorev2.get_icd_codes`` and ``pyPheWAS.pyPhewasCorev2.get_cpt_codes``. If a subject in
           the group *does not have any records in the phenotype DataFrame* (e.g. if none of their ICD/CPT codes had a
           correponding phenotype), they are *still included* in the feature matrix.

	"""

	# sort genotype and phenotype dataframes by 'id'
	print('Sorting phenotype data...')
	phenotype.sort_values(by=['id','AgeAt'+code_type], inplace=True)
	print('Sorting group data...')
	genotypes.sort_values(by='id', inplace=True)
	genotypes.reset_index(inplace=True,drop=True)

	# use phewas/prowas codes to make a dictionary of indices in the np array
	phecode_col = pheno_map[code_type]['map_key']
	empty_pheno_df = pheno_map[code_type]['codes'].set_index(phecode_col)

	empty_pheno_df.sort_index(inplace=True)
	empty_pheno_df['np_index'] = range(0,empty_pheno_df.shape[0])
	np_index = empty_pheno_df['np_index'].to_dict()

	feature_matrix = np.zeros((3, genotypes.shape[0], empty_pheno_df.shape[0]), dtype=float)

	# make genotype a dictionary for faster access time
	genotypes_dict = genotypes.set_index('id').to_dict('index')

	exclude = []  # list of ids to exclude (in icd list but not in genotype list)
	last_id = ''  # track last id seen in icd list
	count = -1
	age_col = 'MaxAgeAt'+code_type

	for _,event in tqdm(phenotype.iterrows(), desc="Processing "+code_type+"s", total=phenotype.shape[0]):
		curr_id = event['id']
		if not curr_id in genotypes_dict:
			if not curr_id in exclude:
				print('%s has records in phenotype file but is not in group file - excluding from study' % curr_id)
				exclude.append(curr_id)
			continue
		# check id to see if a new subject has been found
		if last_id != curr_id:
			count += 1
			while(curr_id != genotypes.loc[count,'id']): # subject at genotypes.loc[count] does not have ICD codes
				empty_id = genotypes.loc[count,'id']
				feature_matrix[1][count] = genotypes_dict[empty_id]['MaxAgeAtVisit']
				count += 1 # continue processing next subject
			last_id = curr_id  # reset last_id
			feature_matrix[1][count] = genotypes_dict[curr_id]['MaxAgeAtVisit']

		# get column index of event phecode
		phecode_ix = np_index[event[phecode_col]]
		# add the max age of that ICD code to the age feature matrix (events are pre-sorted by age)
		feature_matrix[1][count][phecode_ix] = event[age_col]
		# if using a phecode covariate and the event phecode is equal to the covariate phecode, set the fm2 to 1 for this subject
		if (pheno_cov is not None) & (event[phecode_col] == pheno_cov):
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

	return feature_matrix,list(empty_pheno_df.index)



"""

Statistical Modeling

"""


def get_phenotype_info(p_index, code_type):
	"""
	Retrieve phenotype info.

	Return info for the phenotype code at the given index in the phenotype map. The phenotype map (PheWAS or ProWAS)
	depends on the value of ``code_type``. For ``code_type`` = 'ICD', PheWAS code info is returned. For ``code_type`` =
	'CPT', ProWAS code info is returned.

	The returned information consists of (in order) PheWAS/ProWAS code, Phenotype, and ICD-9/ICD-10 or CPT rollup (i.e. all
	ICD/CPT codes that map to the given PheWAS/ProWAS code)

	:param p_index: index of the desired phenotype code
	:param code_type: type of phenotype map to use
	:type p_index: int
	:type code_type: str

	:returns: list of phenotype info
	:rtype: list of str
	"""
	code_map = pheno_map[code_type]['codes']

	p_code = code_map.loc[p_index, pheno_map[code_type]['map_key']]
	p_name = code_map.loc[p_index, pheno_map[code_type]['map_name']]
	p_cat = code_map.loc[p_index, pheno_map[code_type]['cat_name']]

	if code_type == 'ICD':
		icd9_ix = icd9_codes['PheCode'] == p_code
		icd10_ix = icd10_codes['PheCode'] == p_code
		p_icd9_rollup = '/'.join(icd9_codes.loc[icd9_ix, 'ICD_CODE'].tolist())
		p_icd10_rollup = '/'.join(icd10_codes.loc[icd10_ix, 'ICD_CODE'].tolist())
		return [p_code, p_name, p_cat, p_icd9_rollup, p_icd10_rollup]
	else:
		cpt_ix = cpt_codes['prowas_code'] == p_code
		p_cpt_rollup = '/'.join(cpt_codes.loc[cpt_ix, 'CPT_CODE'].tolist())
		return [p_code, p_name, p_cat, p_cpt_rollup]



def fit_pheno_model(model_str, model_type, model_data, phe_thresh=5):
	"""
	Compute the specified logistic regression. Note: model intercept automatically added.


	:param model_str: a patsy-like regression formula
	:param model_type: type of model [linear (GLM - Gaussian family w/ identity link), log (logistic with regularized max likelihood)]
	:param model_data: data for estimating model; all variables from the model_str must be included as columns.
	:param phe_thresh: *[optional]* threshold for running regression; see note (default: *5*)

	:type model_str: str
	:type model_type: str
	:type model_data: pandas DataFrame
	:type phe_thresh: int

	:returns: (regression_model, note)
	:rtype: tuple (statsmodels model, str)

	.. note:: To prevent false positives & improve statistical power, regressions are only computed for
		phenotypes which present for greater than ``phe_thresh`` subjects in the cohort.

	"""
	
	note = None

	if not np.where(model_data['phe'] > 0)[0].shape[0] > phe_thresh:
		note = f"ERROR < {phe_thresh} records with phecode"
		model = None
	else:
		try:
			with warnings.catch_warnings(record=True) as cx_manager:
				if model_type == "linear":
					model = smf.glm(model_str, model_data).fit(disp=False)
				else:
					model = smf.logit(model_str, model_data).fit_regularized(method='l1', alpha=0.1, disp=0, trim_mode='size', qc_verbose=0)
				# deal with statsmodels warnings
				if len(cx_manager) > 0: note = '|'.join( [str(w.message).replace('\n','. ') for w in cx_manager] )
		except ValueError as ve:
			note = f"ERROR computing regression: {str(ve).strip()}"
			model = None
		except Exception as e:
			note = f"ERROR computing regression: {str(e).strip()}"
			model = None
	return (model, note)



def parse_pheno_model(reg, phe_model, note, phe_info, var):
	"""
	Parse results from fit_pheno_model()

	:param reg: regression dataframe
	:param phe_model: regression model for phenotype var
	:param phe_info: metadata for phenotype var
	:param var: variable for which to collect results

	:returns: None
	"""

	ix = reg.shape[0] # next available index in reg
	if phe_model is not None:
		p = phe_model.pvalues[var]
		beta = phe_model.params[var]
		conf = phe_model.conf_int()
		conf_int = '[%s,%s]' % (conf[0][var], conf[1][var])
		stderr = phe_model.bse[var]
		stat_info = [-math.log10(p), p, beta, conf_int, stderr]
		if phe_model.summary().extra_txt:
			note = phe_model.summary().extra_txt.replace('\n', '. ')
	else:
		stat_info = [np.nan, np.nan, np.nan, None, np.nan]

	reg.loc[ix] = phe_info[0:2] + [note] + stat_info + phe_info[2:]
	return



def run_phewas(fm, demo, code_type, reg_type, covariates='', response='genotype', phe_thresh=5, canonical=True):
	"""
	Run mass phenotype regressions

	Iterate over all PheWAS/ProWAS codes in the feature matrix, running a logistic regression of the form:

	:math:`Pr(phenotype\_aggregate) \sim logit(response + covariates)`

	or the *reverse* form (`canonical=False`):

	:math:`Pr(response) \sim logit(phenotype\_aggregate + covariates)`

	``fm`` is a 3xNxP matrix, where N = number of subjects and P = number of PheWAS/ProWAS Codes; this should only
	be consutrcted by ``pyPheWAS.pyPhewasCorev2.generate_feature_matrix`` - otherwise results will be untrustworthy.
	To use the age feature matrix (``fm[1]``), include 'MaxAgeAtICD' or 'MaxAgeAtCPT' in the ``covariates`` string.
	Other than 'MaxAgeAtICD' and 'MaxAgeAtCPT', all covariates and the response variable must be included in
	the group DataFrame.

	The returned DataFrame includes the PheWAS/ProWAS code, Phenotype (code description, e.g. 'Pain in joint'),
	-log\ :sub:`10`\ (p-value), p-value, beta, beta's confidence interval, beta's standard error, and lists of the ICD-9/ICD-10 or
	CPT codes that map to the phenotype.

	:param fm: phenotype feature matrix derived via ``pyPheWAS.pyPhewasCorev2.generate_feature_matrix``
	:param demo: group data
	:param code_type:  type of EMR code ('ICD' or 'CPT')
	:param reg_type: type of regression (0:binary, 1:count, 2:duration)
	:param covariates: *[optional]* covariates to include in the regressions separated by '+' (e.g. 'sex+ageAtDx')
	:param response: *[optional]* response variable in the logisitc model (default: *genotype*)
	:param phe_thresh: *[optional]* threshold for running regression; see note (default: *5*)
	:param canonical: *[optional]*  if False, use the reverse regression formula. if True [default] use the canonical formula.

	:type fm: numpy array
	:type demo: pandas DataFrame
	:type code_type: str
	:type reg_type: int
	:type covariates: str
	:type response: str
	:type phe_thresh: int
	:type reverse: bool

	:returns: tuple -> (regression results, model equation)
	:rtype: tuple(pandas DataFrame, str)


	.. note:: To prevent false positives & improve statistical power, regressions are only computed for
		  phenotypes which present for greater than ``phe_thresh`` subjects in the cohort.

	"""

	# sort group data by id
	print('Sorting group data...')
	demo.sort_values(by='id', inplace=True)
	demo.reset_index(inplace=True, drop=True)


	fm_shape = len(fm[0, 0])
	# check number of phenotypes in the feature matrix
	num_pheno = pheno_map[code_type]['codes'].shape[0]
	assert fm_shape == num_pheno, "Expected %d columns in feature matrix, but found %d. Please check the feature matrix" % (num_pheno, fm_shape)
	
	### define model ###############################################################
	cols = covariates.split('+') + [response]

	age_col = 'MaxAgeAtICD' if (code_type == 'ICD') else 'MaxAgeAtCPT'
	if age_col in cols:
		cols.remove(age_col)
	else:
		age_col = None
	
	model_data = demo[cols].copy()

	phe_cov = 'phewas_cov' if (code_type == 'ICD') else 'prowas_cov'
	if fm[2].any(): # phewas_cov feature matrix
		model_data[phe_cov] = fm[2][:,0] # all columns are the same in this one
		covariates += f'+{phe_cov}'

	if not canonical:
		model_str = f"{response}~phe+{covariates}"
	else:
		model_str = f"phe~{response}+{covariates}"

	model_type = "linear" if canonical and (reg_type != 0) else "log"
	res_var = response if canonical else "phe"
	################################################################################

	# look for perfect separation between response variable and phecodes
	r_0 = fm[0][model_data[response] == 0, :]
	r_1 = fm[0][model_data[response] == 1, :]
	sep_flat = (r_0.any(axis=0) & ~r_1.any(axis=0)) | (~r_0.any(axis=0) & r_1.any(axis=0))

	regressions = pd.DataFrame(columns=pheno_map[code_type]['reg_cols'])

	for index in tqdm(range(fm_shape), desc='Running Regressions'):
		phen_info = get_phenotype_info(index, code_type)

		model_data['phe'] = fm[0][:, index] # aggregate phenotype data 
		if age_col is not None: model_data[age_col] = fm[1][:, index] # MaxAgeAtEvent

		phe_model, note = fit_pheno_model(model_str, model_type, model_data, phe_thresh)
		if sep_flat[index] and (note is None): note = "WARNING perfect separation between phecode and response"

		parse_pheno_model(regressions, phe_model, note, phen_info, res_var)

		model_data.drop(columns=['phe'], inplace=True)
		if age_col is not None: model_data.drop(columns=[age_col], inplace=True)

	return regressions.sort_values(by='p-val'), model_str # sort by significance # TODO: make sure this works!!!! (bin, lin, and dur)


"""

Result Visualization

"""


def get_bon_thresh(p_values, alpha=0.05):
	"""
	Calculate the Bonferroni multiple comparisons correction threshold for a list of p-values.

	Divide the power by the number of finite p-values values (all non-nan values).

	:param p_values: list of p-values
	:param alpha: the uncorrected significance level being used (default = 0.05)
	:type p_values: numpy array
	:type alpha: float

	:returns: The Bonferroni correction threshold
	:rtype: float

	"""
	return alpha / sum(np.isfinite(p_values))


def get_fdr_thresh(p_values, alpha=0.05):
	"""
	Calculate the false discovery rate (FDR) multiple comparisons correction threshold for a list of p-values.

	:param p_values: list of p-values
	:param alpha: the uncorrected significance level being used (default = 0.05)
	:type p_values: numpy array
	:type alpha: float

	:returns: The FDR correction threshold
	:rtype: float

	"""
	sn = np.sort(p_values)
	sn = sn[np.isfinite(sn)]
	for i in range(len(sn)):
		p_crit = alpha * float(i+1) / float(len(sn))
		if sn[i] <= p_crit:
			continue
		else:
			break
	return sn[i]



def get_bhy_thresh(p_values, alpha):
	# Deprecated
	"""
	Calculate the false discovery rate threshold.

	:param p_values: a list of p-values obtained by executing the regression
	:param alpha: the thershold power being used (usually 0.05)
	:type p_values: numpy array
	:type alpha: float

	:returns: the false discovery rate
	:rtype: float
	"""
	sn = np.sort(p_values)
	sn = sn[np.isfinite(sn)]
	# sn = sn[::-1]
	for i in range(len(sn)):
		p_crit = alpha * float(i+1) / (8.1 * float(len(sn)))
		if sn[i] <= p_crit:
			continue
		else:
			break
	return sn[i]


def get_imbalances(regressions):
	# Deprecated
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


def get_x_label_positions(categories, lines=True):
	# Deprecated
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


def plot_manhattan(regressions, thresh, code_type='ICD', show_imbalance=True, plot_all_pts=True, save='', save_format=''):
	"""
	Plots significant phenotype data on a Manhattan Plot.

	The significance of each phenotype (represented by :math:`-log_{10}(p)`\ ) is plotted along the
	y-axis, with phenotypes plotted along the x-axis.
	If ``save`` is provided, the plot is saved to a file; otherwise, the plot may be displayed with
	matplotlib.pyplot.show() after this function returns.


	:param regressions: dataframe containing the regression results
	:param thresh: p-value significance threshold
	:param code_type: type of EMR data ('ICD' or 'CPT')
	:param show_imbalance: boolean variable that determines whether or not to show imbalances on the plot (default True)
	:param plot_all_pts: plot all points regardless of significance (default True)
	:param save: the output file to save to (if empty, display the plot)
	:param save_format: format of the save file
	:type regressions: pandas DataFrame
	:type thresh: float
	:type code_type: str
	:type show_imbalance: boolean
	:type save: str
	:type save_format: str

	"""

	# Initialize figure
	fig = plt.figure(1)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# Merge regressions with phenotype data to get categories (if they're not already present)
	if (pheno_map[code_type]['cat_key'] is not None) and not (pheno_map[code_type]['cat_key'] in regressions.columns):
		regressions = pd.merge(regressions, pheno_map[code_type]['codes'],
							   left_on=pheno_map[code_type]['pheno_key'],
							   right_on=pheno_map[code_type]['map_key'],
							   suffixes=[None, '_y'])

	# Sort by category if we can
	if (pheno_map[code_type]['cat_key'] is not None):
		regressions.sort_values(by=pheno_map[code_type]['cat_key'], inplace=True)

	# Plot all points w/ labels
	e = 1 # x-axis counter
	artists = []
	plt.ylabel('-log10(p)')

	ax.axhline(y=-math.log10(thresh), color='red', ls='dotted')  # plot threshold

	for ix,data in regressions.iterrows():
		logp_ix = data['"-log(p)"']
		# determine marker type based on whether/not showing imbalance
		if show_imbalance:
			mew = 1.5
			if data['beta'] > 0: m = '+'
			else: m = '_'
		else:
			mew = 0.0
			m = 'o'
		# get point color
		if (code_type == 'ICD') and (pheno_map[code_type]['cat_key'] is not None):
			c = cat_colors[data[pheno_map[code_type]['cat_name']]]
		else:
			c = 'xkcd:aqua' # constant color
		# Plot PheCode data point & format PheCode label
		if data['p-val'] < thresh:
			# plot significant PheCode/ProCode with label
			ax.plot(e, logp_ix, m, color=c, fillstyle='full', markeredgewidth=mew)
			artists.append(ax.text(e, logp_ix, data[pheno_map[code_type]['pheno_name']],
								   rotation=45, va='bottom', fontsize=6))
			e += 15
		elif plot_all_pts:
			# plot insignificant PheCode/ProCode without label
			ax.plot(e, logp_ix, m, color=c, fillstyle='full', markeredgewidth=mew)
			e += 15
	# Category Legend
	if (code_type == 'ICD') and (pheno_map[code_type]['cat_key'] is not None):
		line1 = []
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
		for lab in cat_colors.keys():
			line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=cat_colors[lab], label=lab))
		artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, 0), loc='upper center', fancybox=True, ncol=4, prop={'size': 6}))

	# Plot x axis
	ax.axhline(y=0, color='black')
	frame1.axes.get_xaxis().set_visible(False)

	# Save the plot
	if save:
		plt.savefig(save, format = save_format, bbox_extra_artists = artists, bbox_inches ='tight', dpi = 300)
		plt.close()
	
	return


def plot_log_odds_ratio(regressions, thresh, code_type='ICD', save='', save_format='', label_loc="plot"):
	"""
	Plots significant phenotype data on a Log Odds Plot.

	The log odds value & confidence interval is plotted along the x-axis, with Phenotypes sorted by category
	plotted along the y-axis.
	If ``save`` is provided, the plot is saved to a file; otherwise, the plot may be displayed with
	matplotlib.pyplot.show() after this function returns.

	:param regressions: dataframe containing the regression results
	:param thresh: p-value significance threshold
	:param code_type: type of EMR data ('ICD' or 'CPT')
	:param save: the output file to save to (if empty, display the plot)
	:param save_format: format of the save file
	:param label_loc: where to plot Phenotype labels ["plot" (defulat) or "axis"]
	:type regressions: pandas DataFrame
	:type thresh: float
	:type code_type: str
	:type save: str
	:type save_format: str
	:type label_loc: str

	"""

	# Initialize figure
	fig = plt.figure(2)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# Merge regressions with phenotype data to get categories (if they're not already present)
	if (pheno_map[code_type]['cat_key'] is not None) and not (pheno_map[code_type]['cat_key'] in regressions.columns):
		regressions = pd.merge(regressions, pheno_map[code_type]['codes'],
							   left_on=pheno_map[code_type]['pheno_key'],
							   right_on=pheno_map[code_type]['map_key'],
							   suffixes=[None, '_y'])

	# Sort by category if we can
	if (pheno_map[code_type]['cat_key'] is not None):
		regressions.sort_values(by=pheno_map[code_type]['cat_key'], inplace=True)

	# Plot all points w/ labels
	e = 1 # vertical index
	text_size = 6
	artists = []
	pheno_name = pheno_map[code_type]['pheno_name']
	if label_loc == "axis":
		phecode_labels = []
		phecode_locs = []
	plt.xlabel('Log odds ratio')
	for ix, data in regressions.iterrows():
		beta_ix = data['beta']
		if  data['p-val'] < thresh:
			# Add Phecode label
			if label_loc == "plot":
				if beta_ix > 0: # only difference is ha (horizontal alignment)
					artists.append(ax.text(beta_ix, e, data[pheno_name], rotation=0, ha='left', fontsize=text_size))
				else:
					artists.append(ax.text(beta_ix, e, data[pheno_name], rotation=0, ha='right', fontsize=text_size))
			else: # location = "axis"
				phecode_labels.append(data[pheno_name])
				phecode_locs.append(e)
			# get point color
			if (code_type == 'ICD') and (pheno_map[code_type]['cat_key'] is not None):
				c = cat_colors[data[pheno_map[code_type]['cat_name']]]
			else:
				c = 'xkcd:aqua' # constant color

			# Plot Phecode Data
			ax.plot(beta_ix, e, 'o', color=c, fillstyle='full', markeredgewidth=0.0)
			ax.plot([data['lowlim'], data['uplim']], [e, e], color=c)
			e += 15

	# Plot y axis
	ax.axvline(x=0, color='black')

	if label_loc == "axis":
		plt.yticks(phecode_locs,phecode_labels, ha='right',fontsize=text_size)
	else:
		frame1.axes.get_yaxis().set_visible(False)

	# Legend
	if (code_type == 'ICD') and (pheno_map[code_type]['cat_key'] is not None):
		line1 = []
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
		for lab in cat_colors.keys():
			line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=cat_colors[lab], label=lab))
		artists.append(ax.legend(handles=line1, bbox_to_anchor=(0.5, -0.125), loc='upper center', fancybox=True, ncol=4, prop={'size': text_size}))


	# Save the plot
	if save:
		plt.savefig(save,
					format = save_format,
					bbox_extra_artists = artists,
					bbox_inches ='tight',
					dpi = 300
					)
		plt.close()

	return


def plot_volcano(regressions, code_type='ICD', save='', save_format=''):
	"""
	Plots all phenotype data on a Volcano Plot.

	The significance of each phenotype (represented by :math:`-log_{10}(p)`\ ) is plotted along the
	y-axis, with log odds value (effect size) plotted along the x-axis. To improve plot legibility, only those
	phenotypes which surpass the FDR/Bonferroni significance thresholds have labels displayed.
	If ``save`` is provided, the plot is saved to a file; otherwise, the plot may be displayed with
	matplotlib.pyplot.show() after this function returns.

	:param regressions: dataframe containing the regression results
	:param code_type: type of EMR data ('ICD' or 'CPT')
	:param save: the output file to save to (if empty, display the plot)
	:param save_format: format of the save file
	:type regressions: pandas DataFrame
	:type code_type: str
	:type save: str
	:type save_format: str

	"""

	# get thresh values
	bon = get_bon_thresh(regressions["p-val"].values, 0.05)
	fdr = get_fdr_thresh(regressions["p-val"].values, 0.05)

	# Initialize figure
	fig = plt.figure(3)
	ax = plt.subplot(111)
	frame1 = plt.gca()

	# Plot all points w/ labels
	artists = []
	plt.ylabel('-log10(p)')
	plt.xlabel('Log Odds Ratio')

	regressions.sort_values(by='p-val', ascending=False, inplace=True)

	for ix,data in regressions.iterrows():
		logp_ix = data['"-log(p)"']
		beta = data['beta']
		# determine marker color & label based on thresholds
		if  data['p-val'] < bon:
			c = 'gold'
			phe = data[pheno_map[code_type]['pheno_name']]
		elif data['p-val'] < fdr:
			c = 'midnightblue'
			phe = data[pheno_map[code_type]['pheno_name']]
		else:
			c = 'slategray'
			phe = ''

		# Plot PheCode data point & format PheCode label
		ax.plot(beta, logp_ix, 'o', color=c, fillstyle='full', markeredgewidth=0)
		artists.append(ax.text(beta, logp_ix, phe, rotation=45, va='bottom', fontsize=4))

	# Legend
	line1 = []
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])

	line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="gold", label="BonFerroni"))
	line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="midnightblue", label="FDR"))
	line1.append(mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="slategray", label="Insignificant"))

	artists.append(ax.legend(handles=line1, loc='best', fancybox=True, ncol=4, prop={'size': 6}))

	# Save the plot
	if save:
		plt.savefig(save,
					format = save_format,
					bbox_extra_artists = artists,
					bbox_inches ='tight',
					dpi = 300
					)
		plt.close()

	return



"""

Housekeeping & Parameters

"""


def process_args(kwargs, optargs, *args):
	# Deprecated
	clean = np.vectorize(lambda x: x[x.rfind('-') + 1:] + '=')
	searchfor = clean(list(optargs.keys()))
	opts, rem = getopt.getopt(args, '', searchfor)
	assert len(rem) == 0, 'Unknown arguments included %s' % (str(rem))
	for option in opts:
		k, v = option
		kwargs[optargs[k]] = v

	return kwargs


def display_kwargs(kwargs):
	for k, v in kwargs.items():
		num_dots = 80 - len(str(k)) - len(str(v))
		print(str(k) + '.'*num_dots + str(v))


cat_colors = {'other': 'gold',
			  'circulatory system': 'xkcd:bright red',
			  'congenital anomalies': 'mediumspringgreen',
			  'dermatologic': 'xkcd:dark peach',
			  'digestive': 'yellowgreen',
			  'endocrine/metabolic': 'darkred',
			  'genitourinary': 'seagreen',
			  'hematopoietic': 'orange',
			  'infectious diseases': 'blue',
			  'injuries & poisonings': 'slategray',
			  'mental disorders': 'xkcd:hot pink',
			  'musculoskeletal': 'darkkhaki',
			  'neoplasms': 'xkcd:bluish',
			  'neurological': 'xkcd:purplish pink',
			  'pregnancy complications': 'peachpuff',
			  'respiratory': 'xkcd:carolina blue',
			  'sense organs': 'darkviolet',
			  'symptoms': 'aqua'}

old_plot_colors = {'-': 'gold',
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

# make plot text editable
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


RESERVED_COL_NAMES = ['MaxAgeAtICD', 'MaxAgeAtCPT', 'phe', 'phewas_cov', 'prowas_cov']

#----------------------------------------------------------
# load ICD maps (pyPheWAS)
icd9_codes = get_codes('phecode_map_v1_2_icd9.csv')
icd10_codes = get_codes('phecode_map_v1_2_icd10_beta.csv')
# load CPT maps (pyProWAS)
cpt_codes = get_codes('prowas_codes.csv')
#----------------------------------------------------------

# Get PheCode list (merge ICD9 & ICD10 maps)
if sum(icd9_codes.columns.isin(['category','category_string'])) == 2:
	mcols = ['PheCode','Phenotype','category','category_string']
	phe_cats = True
else:
	mcols = ['PheCode','Phenotype']
	phe_cats = False
phewas_codes = pd.concat([icd9_codes, icd10_codes[mcols]], sort=False, ignore_index=True)
phewas_codes = phewas_codes[mcols].dropna()
phewas_codes.drop_duplicates(subset='PheCode',inplace=True)
phewas_codes.sort_values(by=['PheCode'], inplace=True)
phewas_codes.reset_index(inplace=True, drop=True)

# Get ProCode list
if sum(cpt_codes.columns.isin(['ccs','CCS Label'])) == 2:
	mcols = ['prowas_code','prowas_desc','ccs','CCS Label']
	pro_cats = True
else:
	mcols = ['prowas_code','prowas_desc']
	pro_cats = False
prowas_codes = cpt_codes[mcols].drop_duplicates(subset='prowas_code')
prowas_codes.sort_values(by=['prowas_code'], inplace=True)
prowas_codes.reset_index(inplace=True,drop=True)


# Build Phenotype Information Map
pheno_map = {}
pheno_map['ICD'] = {'codes' : phewas_codes,
					'pheno_name' : 'PheWAS Name',
					'pheno_key' : 'PheWAS Code',
					'map_key' : 'PheCode',
					'map_name' : 'Phenotype',
					'cat_key' : 'category' if phe_cats else None,
					'cat_name' : 'category_string' if phe_cats else None,
					'reg_cols' : ['PheWAS Code', 'PheWAS Name', 'note', '\"-log(p)\"',
								  'p-val', 'beta', 'Conf-interval beta',
								  'std_error', 'category_string', 'ICD-9', 'ICD-10']
}

pheno_map['CPT'] = {'codes' : prowas_codes,
					'pheno_name' : 'ProWAS Name',
					'pheno_key' : 'ProWAS Code',
					'map_key' : 'prowas_code',
					'map_name' : 'prowas_desc',
					'cat_key' : 'ccs' if pro_cats else None,
					'cat_name' : 'CCS Label' if pro_cats else None,
					'reg_cols' : ['ProWAS Code', 'ProWAS Name', 'note', '\"-log(p)\"',
								  'p-val', 'beta', 'Conf-interval beta',
								  'std_error', 'CCS Label', 'CPT']
}
