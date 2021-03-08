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
	phenotypes = pd.merge(cptfile, cpt_codes, on='CPT_CODE')

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
	phenotype.sort_values(by='id', inplace=True)
	print('Sorting group data...')
	genotypes.sort_values(by='id', inplace=True)
	genotypes.reset_index(inplace=True,drop=True)

	# use phewas/prowas codes to make a dictionary of indices in the np array
	if code_type == 'CPT':
		phecode_col = 'prowas_code'
		empty_phewas_df = prowas_codes.set_index('prowas_code')
	else: # code_type == 'ICD'
		phecode_col = 'PheCode'
		empty_phewas_df = phewas_codes.set_index('PheCode')

	empty_phewas_df.sort_index(inplace=True)
	empty_phewas_df['np_index'] = range(0,empty_phewas_df.shape[0])
	np_index = empty_phewas_df['np_index'].to_dict()

	feature_matrix = np.zeros((3, genotypes.shape[0], empty_phewas_df.shape[0]), dtype=float)

	# make genotype a dictionary for faster access time
	genotypes_dict = genotypes.set_index('id').to_dict('index')

	exclude = []  # list of ids to exclude (in icd list but not in genotype list)
	last_id = ''  # track last id seen in icd list
	count = -1

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
		# add the max at of that ICD code to the age feature matrix
		feature_matrix[1][count][phecode_ix] = event['MaxAgeAt'+code_type]
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

	return feature_matrix,list(empty_phewas_df.index)



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
	if code_type == 'ICD':
		p_code = phewas_codes.loc[p_index, 'PheCode']
		p_name = phewas_codes.loc[p_index, 'Phenotype']
		p_cat = phewas_codes.loc[p_index, 'category_string']
		icd9_ix = icd9_codes['PheCode'] == p_code
		icd10_ix = icd10_codes['PheCode'] == p_code

		p_icd9_rollup = '/'.join(icd9_codes.loc[icd9_ix, 'ICD_CODE'].tolist())
		p_icd10_rollup = '/'.join(icd10_codes.loc[icd10_ix, 'ICD_CODE'].tolist())

		return [p_code, p_name, p_cat, p_icd9_rollup, p_icd10_rollup]

	else:
		p_code = prowas_codes.loc[p_index, 'prowas_code']
		p_name = prowas_codes.loc[p_index, 'prowas_desc']
		p_cat = prowas_codes.loc[p_index, 'CCS Label']
		cpt_ix = cpt_codes['prowas_code'] == p_code
		p_cpt_rollup = '/'.join(cpt_codes.loc[cpt_ix, 'CPT_CODE'].tolist())

		return [p_code, p_name, p_cat, p_cpt_rollup]


def fit_pheno_model(genotypes, phen_vector1, phen_vector2, phen_vector3='', covariates='',
                         response='genotype',  lr=0, phenotype='', code_type='ICD'):
	"""
	Runs a logistic regression for a specific phenotype vector

	Compute a logistic regression of the form:

	:math:`Pr(response) \sim logit(phen\_vector1 + covariates)`

	``phen_vector1`` is a vector of phenotype aggregates of length N, where N = number of subjects.
	To use the age feature vector (``phen_vector2``), include 'MaxAgeAtICD' or 'MaxAgeAtCPT' in the ``covariates`` string.
	Other than 'MaxAgeAtICD' and 'MaxAgeAtCPT', all covariates and the response variable must be included in
	the group DataFrame.

	The returned results list consists of (in order) the -log\ :sub:`10`\ (p-value), p-value, beta, beta's confidence interval,
	and beta's standard error, estimated from the logit model for the ``phen_vector1`` variable.

	:param genotypes: group data
	:param phen_vector1: the aggregate phenotype vector
	:param phen_vector2: the maximum event age phenotype vector
	:param phen_vector3: *[optional]* the phenotype covariate vector
	:param covariates: *[optional]* covariates to include in the regressions separated by '+' (e.g. 'sex+ageAtDx')
	:param response: *[optional]* response variable in the logit model (default: *genotype*)
	:param lr: *[optional]* regularized maximum likelihood optimization flag (0 [default] = not regularized; 1 = regularized)
	:param phenotype: *[optional]* phenotype info [code, description] for this regression (used only for error handling)
	:param code_type: *[optional]*  EMR data type ('ICD' [default], or 'CPT')
	:type genotypes: pandas DataFrame
	:type phen_vector1: numpy array
	:type phen_vector2: numpy array
	:type phen_vector3: numpy array
	:type covariates: str
	:type response: str
	:type lr: int [0,1]
	:type phenotype: list of str
	:type code_type: str

	:returns: regression results
	:rtype: list

	"""

	data = genotypes.copy()
	data['y'] = phen_vector1 # aggregate phenotype data

	# append '+' to covariates (if there are any) -> makes definition of 'f' more elegant
	if covariates != '':
		covariates = '+' + covariates

	# add MaxAgeAtEvent to data (if needed)
	if (code_type == 'ICD') and ('MaxAgeAtICD' in covariates):
		data['MaxAgeAtICD'] = phen_vector2
	elif (code_type == 'CPT') and ('MaxAgeAtCPT' in covariates):
		data['MaxAgeAtCPT'] = phen_vector2

	# add phewas_cov feature matrix to data & covariates (if needed)
	if phen_vector3.any():
		data['phe_cov'] = phen_vector3
		covariates = covariates + '+ phe_cov'

	# define model ('f') for the logisitc regression
	if lr == 1:
		predictors = covariates.replace(" ", "").split('+')
		predictors[0] = 'y'
		f = [response.strip(), predictors]
	else:
		f = response + '~ y' + covariates

	try:
		if lr == 0:
			# fit logit without regulatization
			model = smf.logit(f, data).fit(disp=False)
		elif lr == 1:
			# fit logit with regularization
			logit = sm.Logit(data[f[0]], data[f[1]])
			model = logit.fit_regularized(method='l1', alpha=0.1, disp=0, trim_mode='size', qc_verbose=0)
		else:
			linreg = smf.logit(f, data).fit(method='bfgs', disp=False)
		# get results
		p = model.pvalues.y
		beta = model.params.y
		conf = model.conf_int()
		conf_int = '[%s,%s]' % (conf[0]['y'], conf[1]['y'])
		stderr = model.bse.y
		reg_result = [-math.log10(p), p, beta, conf_int, stderr]  # collect results

	except ValueError as ve:
		print('\n')
		if phenotype != '':
			print('ERROR computing regression for phenotype %s (%s)' %(phenotype[0],phenotype[1]))
		print(ve)
		print('lr = % d' %lr)
		reg_result = [np.nan, np.nan, np.nan, np.nan, np.nan]
	except Exception as e:
		print('\n')
		if phenotype != '':
			print('ERROR computing regression for phenotype %s (%s)' %(phenotype[0],phenotype[1]))
		print(e)
		reg_result = [np.nan, np.nan, np.nan, np.nan, np.nan]
	return reg_result


def run_phewas_legacy(fm, genotypes, code_type, covariates='', response='genotype', phe_thresh=5):
	"""
	Run mass phenotype regressions (legacy version)

	Iterate over all PheWAS/ProWAS codes in the feature matrix, running a logistic regression of the form:

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
	:param genotypes: group data
	:param code_type:  type of EMR code ('ICD' or 'CPT')
	:param covariates: *[optional]* covariates to include in the regressions separated by '+' (e.g. 'sex+ageAtDx')
	:param response: *[optional]* response variable in the logisitc model (default: *genotype*)
	:param phe_thresh: *[optional]* threshold for running regression; see note (default: *5*)
	:type fm: numpy array
	:type genotypes: pandas DataFrame
	:type code_type: str
	:type covariates: str
	:type response: str
	:type phe_thresh: int

	:returns: regression results for each phenotype
	:rtype: pandas DataFrame


	.. note:: To prevent false positives & improve statistical power, regressions are only computed for
		  phenotypes which present for greater than ``phe_thresh`` subjects in the cohort. Phenotypes which do not meet
		  this criteria are not included in the returned regression results.

	.. note:: For phenotypes that present in both the case (``response`` = 1) and control (``response`` = 0) groups, maximum
		  likelihood optimization is used to compute the logistic regression. For phenotypes that only present in
		  one of those groups, regularized maximum likelihood optimization is used.

	"""

	# sort group data by id
	print('Sorting group data...')
	genotypes.sort_values(by='id', inplace=True)
	genotypes.reset_index(inplace=True,drop=True)

	num_pheno = len(fm[0, 0])
	# check number of phenotypes in the feature matrix
	if code_type == 'ICD':
		assert num_pheno == phewas_codes.shape[0], "Expected %d columns in PheWAS feature matrix, but found %d. Please check the feature matrix" %(phewas_codes.shape[0], num_pheno)
	else:
		assert num_pheno == prowas_codes.shape[0], "Expected %d columns in ProWAS feature matrix, but found %d. Please check the feature matrix" %(prowas_codes.shape[0], num_pheno)

	# store all of the pertinent data from the regressions
	if code_type == 'ICD':
		out_cols = phewas_reg_cols
	else:
		out_cols = prowas_reg_cols
	regressions = pd.DataFrame(columns=out_cols)

	control = fm[0][genotypes[response] == 0, :]
	disease = fm[0][genotypes[response] == 1, :]
	# find all phecodes that only present for a single group (ie only controls or only case show the phecode) -> have to use regularization
	inds = np.where((control.any(axis=0) & ~disease.any(axis=0)) | (~control.any(axis=0) & disease.any(axis=0)))[0]
	for index in tqdm(range(num_pheno), desc='Running Regressions'):
		phen_info = get_phenotype_info(index, code_type)
		phen_vector1 = fm[0][:, index]
		phen_vector2 = fm[1][:, index]
		phen_vector3 = fm[2][:, index]
		# to prevent false positives, only run regressions if more than phe_thresh records have positive values
		if np.where(phen_vector1 > 0)[0].shape[0] > phe_thresh:
			if index in inds: # decide whether or not to use regularization in logistic regression
				use_regular = 1
			else:
				use_regular = 0
			stat_info = fit_pheno_model(genotypes, phen_vector1, phen_vector2, phen_vector3, covariates,
			                                 response, phenotype = phen_info[0:2], lr=use_regular, code_type=code_type)
		else:
			# not enough samples to run regression
			stat_info = [np.nan, np.nan, np.nan, np.nan, np.nan]

		# save regression data
		info = phen_info[0:2] + stat_info + phen_info[2:]

		regressions.loc[index] = info

	return regressions.dropna(subset=['p-val']).sort_values(by='p-val') # sort by significance



def run_phewas(fm, genotypes, code_type, covariates='', response='genotype', phe_thresh=5):
	"""
	Run mass phenotype regressions

	Iterate over all PheWAS/ProWAS codes in the feature matrix, running a logistic regression of the form:

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
	:param genotypes: group data
	:param code_type:  type of EMR code ('ICD' or 'CPT')
	:param covariates: *[optional]* covariates to include in the regressions separated by '+' (e.g. 'sex+ageAtDx')
	:param response: *[optional]* response variable in the logisitc model (default: *genotype*)
	:param phe_thresh: *[optional]* threshold for running regression; see note (default: *5*)
	:type fm: numpy array
	:type genotypes: pandas DataFrame
	:type code_type: str
	:type covariates: str
	:type response: str
	:type phe_thresh: int

	:returns: regression results for each phenotype
	:rtype: pandas DataFrame


	.. note:: To prevent false positives & improve statistical power, regressions are only computed for
		  phenotypes which present for greater than ``phe_thresh`` subjects in the cohort. Phenotypes which do not meet
		  this criteria are not included in the returned regression results.

	"""

	# sort group data by id
	print('Sorting group data...')
	genotypes.sort_values(by='id', inplace=True)
	genotypes.reset_index(inplace=True, drop=True)

	num_pheno = len(fm[0, 0])
	# check number of phenotypes in the feature matrix
	if code_type == 'ICD':
		assert num_pheno == phewas_codes.shape[0], "Expected %d columns in PheWAS feature matrix, but found %d. Please check the feature matrix" % (phewas_codes.shape[0], num_pheno)
	else:
		assert num_pheno == prowas_codes.shape[0], "Expected %d columns in ProWAS feature matrix, but found %d. Please check the feature matrix" % (prowas_codes.shape[0], num_pheno)

	# store all of the pertinent data from the regressions
	if code_type == 'ICD':
		out_cols = phewas_reg_cols
	else:
		out_cols = prowas_reg_cols
	regressions = pd.DataFrame(columns=out_cols)

	for index in tqdm(range(num_pheno), desc='Running Regressions'):
		phen_info = get_phenotype_info(index, code_type)
		phen_vector1 = fm[0][:, index]
		phen_vector2 = fm[1][:, index]
		phen_vector3 = fm[2][:, index]
		# to prevent false positives, only run regressions if more than phe_thresh records have positive values
		if np.where(phen_vector1 > 0)[0].shape[0] > phe_thresh:
			stat_info = fit_pheno_model(genotypes, phen_vector1, phen_vector2, phen_vector3, covariates,
										response, phenotype=phen_info[0:2], lr=1, code_type=code_type)
		else:
			# not enough samples to run regression
			stat_info = [np.nan, np.nan, np.nan, np.nan, np.nan]

		# save regression data
		info = phen_info[0:2] + stat_info + phen_info[2:]

		regressions.loc[index] = info

	return regressions.dropna(subset=['p-val']).sort_values(by='p-val')  # sort by significance


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

	if code_type == 'ICD':
		pheno_codes = phewas_codes
		reg_key = 'PheWAS Code'
		map_key = 'PheCode'
		cat_key = 'category'
		cat_label = 'category_string'
		pheno_name = 'PheWAS Name'
	else:
		pheno_codes = prowas_codes
		reg_key = 'ProWAS Code'
		map_key = 'prowas_code'
		cat_key = 'ccs'
		cat_label = 'CCS Label'
		pheno_name = 'ProWAS Name'

	if not cat_key in regressions.columns:
		# Merge regressions with phenotype data to get categories
		regressions = pd.merge(regressions, pheno_codes, left_on=reg_key, right_on=map_key, suffixes=[None, '_y'])

	regressions.sort_values(by=cat_key, inplace=True)

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
		# Plot PheCode data point & format PheCode label
		if code_type == 'ICD': c = cat_colors[data[cat_label]]
		else: c = 'xkcd:aqua' # constant color
		if data['p-val'] < thresh:
			# plot significant PheCode/ProCode with label
			ax.plot(e, logp_ix, m, color=c, fillstyle='full', markeredgewidth=mew)
			artists.append(ax.text(e, logp_ix, data[pheno_name], rotation=45, va='bottom', fontsize=6))
			e += 15
		elif plot_all_pts:
			# plot insignificant PheCode/ProCode without label
			ax.plot(e, logp_ix, m, color=c, fillstyle='full', markeredgewidth=mew)
			e += 15

	# Category Legend
	if code_type == 'ICD':
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
		plt.savefig(save,
					format = save_format,
					bbox_extra_artists = artists,
					bbox_inches ='tight',
					dpi = 300
					)
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

	if code_type == 'ICD':
		pheno_codes = phewas_codes
		reg_key = 'PheWAS Code'
		map_key = 'PheCode'
		cat_key = 'category'
		cat_label = 'category_string'
		pheno_name = 'PheWAS Name'
	else:
		pheno_codes = prowas_codes
		reg_key = 'ProWAS Code'
		map_key = 'prowas_code'
		cat_key = 'ccs'
		cat_label = 'CCS Label'
		pheno_name = 'ProWAS Name'

	if not cat_key in regressions.columns:
		# Merge regressions with phenotype data to get categories
		regressions = pd.merge(regressions, pheno_codes, left_on=reg_key, right_on=map_key, suffixes=[None, '_y'])

	regressions.sort_values(by=cat_key, inplace=True)


	# Plot all points w/ labels
	e = 1 # vertical index
	text_size = 6
	artists = []
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

			# Plot Phecode Data
			if code_type == 'ICD': c = cat_colors[data[cat_label]]
			else: c = 'xkcd:aqua' # constant color
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
	if code_type == 'ICD':
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

	if code_type == 'ICD':
		pheno_name = 'PheWAS Name'
	else:
		pheno_name = 'ProWAS Name'

	for ix,data in regressions.iterrows():
		logp_ix = data['"-log(p)"']
		beta = data['beta']
		# determine marker color & label based on thresholds
		if  data['p-val'] < bon:
			c = 'gold'
			phe = data[pheno_name]
		elif data['p-val'] < fdr:
			c = 'midnightblue'
			phe = data[pheno_name]
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


phewas_reg_cols = ['PheWAS Code',
				  'PheWAS Name',
				  '\"-log(p)\"',
				  'p-val',
				  'beta',
				  'Conf-interval beta',
				  'std_error',
				  'category_string',
				  'ICD-9',
				  'ICD-10']

prowas_reg_cols = ['ProWAS Code',
				  'ProWAS Name',
				  '\"-log(p)\"',
				  'p-val',
				  'beta',
				  'Conf-interval beta',
				  'std_error',
				  'CCS Label',
				  'CPT']

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


# load ICD maps (pyPheWAS)
icd9_codes = get_codes('phecode_map_v1_2_icd9.csv')
# icd9_codes = get_codes('phecode_map_v1_1_icd9.csv') # original ICD-PheCode mapping
icd10_codes = get_codes('phecode_map_v1_2_icd10_beta.csv')

phewas_codes = icd9_codes.append(icd10_codes[['PheCode','Phenotype','category','category_string']])
phewas_codes = phewas_codes[['PheCode','Phenotype','category','category_string']].dropna()
phewas_codes.drop_duplicates(subset='PheCode',inplace=True)
phewas_codes.sort_values(by=['PheCode'], inplace=True)
phewas_codes.reset_index(inplace=True, drop=True)

# load CPT maps (pyProWAS)
cpt_codes = get_codes('prowas_codes.csv')
prowas_codes = cpt_codes[['prowas_code','prowas_desc','ccs','CCS Label']].drop_duplicates(subset='prowas_code')
prowas_codes.sort_values(by=['prowas_code'], inplace=True)
prowas_codes.reset_index(inplace=True,drop=True)
