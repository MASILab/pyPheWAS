"""
PyPheWAS Explorer Core Functions
Developed by:
    Cailey Kerley

MASI Lab
Department of Electrical Engineering and Computer Science
Vanderbilt University
"""

from pyPheWAS.pyPhewasCorev2 import icd9_codes, icd10_codes, phewas_codes, fit_pheno_model
import math
import numpy as np
import os
import pandas as pd
import statsmodels.formula.api as smf
from tqdm import tqdm
import scipy.stats
import http.server
import socketserver
from pathlib import Path


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



def process_group_vars(df):		
	gvars = df.columns.drop('id')
	# fix column names and drop non-numeric columns
	drop_list = []
	gvars_fixed = {}
	for g in gvars:
		if df.dtypes[g] in [int, float, bool]:
			g_new = g[:15] if (len(g) > 15) else g
			g_new = g_new.replace(' ', '_')
			if g_new != g:
				gvars_fixed[g] = g_new
		else:
			drop_list.append(g)
			print(f"WARNING: data type {df.dtypes[g]} in column {g} not supported")
	if len(gvars_fixed.keys()) > 0:
		df.rename(columns=gvars_fixed, inplace=True)
	if len(drop_list) > 0:
		df.drop(columns=drop_list, inplace=True)
	
	gvars = df.columns.drop('id')
	if len(gvars) > 10: 
		print('WARNING: more than 10 group variables found. Limiting display to first 10 group variables')
		df.drop(columns=gvars[10:], inplace=True)
		gvars = df.columns.drop('id')
	
	# Find potential target variables
	drop_list = []
	target_vals = {var : list(uv) for var in gvars if (df.dtypes[var] in [int, bool]) & (len(uv := df[var].unique()) == 2)}
	for var, values in target_vals.items(): # make sure values are 0/1
		if df.dtypes[var] == bool:
			df[var] = df[var].astype(int)
		else:
			check01 = set(values + [0,1])
			if len(check01) == 3:
				fix = [0,1] if (0 in values) else [1,0]
				df.loc[df[var]!= fix[0], var] = fix[1]
			elif len(check01) == 4:
				one_val = np.max(values)
				var_old = f"{var}_old"
				df.rename(columns={var:var_old}, inplace=True)
				df[var] = 0
				df.loc[df[var_old] == one_val, var] = 1
				drop_list.append(var_old)
	target_vars = list(target_vals.keys())
	assert len(target_vars) > 0, "No valid binary target variables found"
	if len(drop_list) > 0:
		df.drop(columns=drop_list, inplace=True)
		gvars = df.columns.drop('id')

	return target_vars, gvars



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
		raise Exception('Found an ICD_TYPE that was not 9 or 10 - Please check phenotype file.')

	# merge with Phecode table, depends on which type(s) of ICD codes are in the icd file
	if icd_types.shape[0] == 2:
		print('Found both ICD-9 and ICD-10 codes.')
		icd9s = icdfile[icdfile['ICD_TYPE'] == 9]
		icd10s = icdfile[icdfile['ICD_TYPE'] == 10]
		phenotypes_9 = pd.merge(icd9s, icd9_codes,on='ICD_CODE')
		phenotypes_10 = pd.merge(icd10s, icd10_codes, on='ICD_CODE')
		phenotypes = pd.concat([phenotypes_9, phenotypes_10], sort=False, ignore_index=True)
	elif icd_types[0] == 9:
		print('Found only ICD-9 codes.')
		phenotypes = pd.merge(icdfile,icd9_codes,on='ICD_CODE')
	elif icd_types[0] == 10:
		print('Found only ICD-10 codes.')
		phenotypes = pd.merge(icdfile, icd10_codes, on='ICD_CODE')
	else:
		raise Exception('An issue occurred while parsing the ICD_TYPE column - Please check phenotype file.')

	# pre-calculate durations for feature matrix step
	phenotypes['duration'] =    phenotypes.groupby(['id', 'PheCode'])['AgeAtICD'].transform('max') - \
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

def get_1D_histogram(group, mask0, var_name):

	bin_max = max(group[var_name])
	bin_min = min(group[var_name])

	H0, edges0 = np.histogram(group.loc[mask0, var_name].to_numpy(),range=(bin_min,bin_max))
	H1, edges1 = np.histogram(group.loc[~mask0, var_name].to_numpy(),range=(bin_min,bin_max))

	H_df0 = pd.DataFrame(index=range(len(H0)))
	H_df0['count'] = H0
	H_df0['xmin'] = edges0[0:-1]
	H_df0['xmax'] = edges0[1:]
	H_df0['response'] = 0

	H_df1 = pd.DataFrame(index=range(len(H1)))
	H_df1['count'] = H1
	H_df1['xmin'] = edges1[0:-1]
	H_df1['xmax'] = edges1[1:]
	H_df1['response'] = 1

	H_df = pd.concat([H_df0, H_df1], sort=False, ignore_index=True)

	H_df['var_name'] = var_name
	return H_df


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
	ix = 0
	# correlation
	res_df.loc[ix, 'test_name'] = 'correlation'
	res_df.loc[ix, 'var'] = 'joint'
	var_data = g_df[[var1, var2]].to_numpy()
	[corr, pval] = scipy.stats.spearmanr(var_data)
	res_df.loc[ix, 'result'] = corr
	res_df.loc[ix, 'pval'] = pval
	ix += 1

	# univariate regressions
	res_df.loc[ix, 'test_name'] = 'univariate regression'
	res_df.loc[ix, 'var'] = var1
	f_v1 = response + ' ~ ' + var1
	logreg1 = smf.logit(f_v1, g_df).fit(disp=False)
	res_df.loc[ix, 'pval'] = logreg1.pvalues[var1]
	res_df.loc[ix, 'result'] = logreg1.params[var1]
	ix+=1
	res_df.loc[ix, 'test_name'] = 'univariate regression'
	res_df.loc[ix, 'var'] = var2
	f_v2 = response + ' ~ ' + var2
	logreg2 = smf.logit(f_v2, g_df).fit(disp=False)
	res_df.loc[ix, 'pval'] = logreg2.pvalues[var2]
	res_df.loc[ix, 'result'] = logreg2.params[var2]
	ix += 1

	# multivariate regression
	f_mul = response + ' ~ ' + var1 + ' + ' + var2
	logreg3 = smf.logit(f_mul, g_df).fit(disp=False)
	res_df.loc[ix, 'test_name'] = 'multivariate regression'
	res_df.loc[ix, 'var'] = var1
	res_df.loc[ix, 'pval'] = logreg3.pvalues[var1]
	res_df.loc[ix, 'result'] = logreg3.params[var1]
	ix += 1
	res_df.loc[ix, 'test_name'] = 'multivariate regression'
	res_df.loc[ix, 'var'] = var2
	res_df.loc[ix, 'pval'] = logreg3.pvalues[var2]
	res_df.loc[ix, 'result'] = logreg3.params[var2]

	return res_df



"""

Statistical Modeling

"""


def get_phenotype_info(p_index): 
	"""
	Returns all of the info of the phewas code at the given index.

	:param p_index: The index of the desired phewas code
	:type p_index: int

	:returns: A list including the code, the name, and the category.
	:rtype: list of strings
	"""
	p_code = phewas_codes.loc[p_index, 'PheCode']
	p_name = phewas_codes.loc[p_index, 'Phenotype']
	p_id = p_name.replace('(','')
	p_id = p_id.replace(')', '')
	p_id = p_id.replace('[', '')
	p_id = p_id.replace(']', '')
	p_id = p_id.replace(',', '')
	p_id = p_id.replace(' ', '')
	p_id = p_id.replace('\'', '')
	p_id = p_id.replace('\"', '')
	p_id = p_id.replace('\\', '')
	p_id = p_id.replace('/', '')

	cat_id = phewas_codes.loc[p_index, 'category']
	cat = phewas_codes.loc[p_index, 'category_string']

	return [p_code, p_name, p_id, cat_id, cat]


def parse_pheno_model(reg, phe_model, phe_info, n_subs, note, var_list):
	"""
	Parse results from fit_pheno_model()

	:param reg: regression dataframe
	:param phe_model: logistic regression model for a phecode 
	:param phe_info: metadata for phecode
	:param n_subs: number of subjects with data present for phecode 
	:param note: notes from model fitting
	:param var_list: list of model variables to save

	:returns: None
	"""

	ix = reg.shape[0] # next available index in reg
	for var in var_list:
		if phe_model is not None:
			p = phe_model.pvalues[var]
			try: log_p = -math.log10(p)
			except:
				log_p = np.nan
				note += "[WARNING p==0.0; PheCode will be excluded from result plots]"
			beta = phe_model.params[var]
			conf = phe_model.conf_int()
			conf_int = '[%s,%s]' % (conf[0][var], conf[1][var])
			stderr = phe_model.bse[var]
			stat_info = [log_p, p, beta, conf_int, conf[0][var], conf[1][var], stderr]

			if phe_model.summary().extra_txt:
				note = "[" + phe_model.summary().extra_txt.replace('\n', '. ') + "]"
		else:
			stat_info = [np.nan, np.nan, np.nan, None, np.nan, np.nan, np.nan]

		reg.loc[ix] = [var] + phe_info[0:3] + [note, n_subs] + stat_info + phe_info[3:]
		ix+=1

	return



def save_pheno_model(reg, var_list, base_path, base_header):
	"""
	Save regression results

	:param reg: regression dataframe
	:param var_list: list of model variables to save
	:param base_path: generic file path for all variables
	:param base_header: generic header for all variables

	:returns: None
	"""
	for var in var_list:
		var_data = reg[reg['result_type'] == var].drop(columns=['result_type'])
		var_file = base_path + var + ("_primary.csv" if var == var_list[0] else '.csv')
		var_header = f"result_variable,{var},{base_header}\n"
		# write regressions to file
		with open(var_file, 'w+') as f:
			f.write(var_header)
			var_data.to_csv(f, index=False)
	return


def check_existing_models(base_path, var_list):
	# check if reg exists - can just reload it instead of re-calculating
	reg_exists = True
	for var in var_list:
		var_file = base_path + var + ("_primary.csv" if var == var_list[0] else '.csv')
		reg_exists &= var_file.exists()
	return reg_exists



def run_phewas(fm, demo, model_str, reg_type, save_cov=False, outpath=Path('.')):
	"""
	For each phewas code in the feature matrix, run the specified type of regression and save all of the resulting stats.

	:param fm: The phewas feature matrix.
	:param demo: A pandas DataFrame containing model covariate and target variables
	:param model_str: a patsy-like regression formula
	:param reg_type: type of regression (0:binary, 1:count, 2:duration)
	:param save_cov: whether/not to save the model's covariate variable results
	:param outpath: path to directory where regression results will be saved

	:returns: A dataframe containing regression data.
	"""

	dependent, predictors = model_str.split('~')
	model_type = "linear" if (dependent == "PheCode") and (reg_type != "binary") else "log"

	base_path = str(outpath/('regressions-%s-%s-%s_' %(reg_type, dependent, predictors)))
	if save_cov:
		var_list = predictors.split('+')
	else:
		var_list = [ predictors.split('+')[0] ]
	
	reg_exists = check_existing_models(base_path, var_list)

	if reg_exists:
		phecode_file = base_path + var_list[0] + '.csv'
		regressions = pd.read_csv(phecode_file, skiprows=1)
	else:
		cols = [dependent] + predictors.split('+')
		cols.remove("PheCode")
		model_data = demo[cols].copy()

		regressions = pd.DataFrame(columns=output_columns)
		num_phecodes = fm.shape[1]

		for index in tqdm(range(num_phecodes), desc='Running Regressions'):
			phe_info = get_phenotype_info(index)
			model_data["PheCode"] = fm[:, index]
			num_nonzero = (model_data["PheCode"] > 0).sum()
			phe_model, note = fit_pheno_model(model_str, model_type, model_data, phe_key = "PheCode")
			parse_pheno_model(regressions, phe_model, phe_info, num_nonzero, note, var_list)
			model_data["PheCode"] = np.nan
		
		if model_type == "log": # Compute Odds Ratio
			regressions['OR'] = np.exp(regressions['beta'])
			regressions['OR_ci_low'] = np.exp(regressions['beta_ci_low'])
			regressions['OR_ci_up'] = np.exp(regressions['beta_ci_up'])

		# Export
		header = ','.join(['model_equation', model_str, 'group', 'group.csv', 'reg_type', reg_type, 'code_type', 'ICD'])
		save_pheno_model(regressions, var_list, base_path, header)
		regressions = regressions[regressions['result_type'] == var_list[0]].drop(columns=['result_type'])

	regressions['pval_str'] = regressions.apply(lambda x: str(x['pval']), axis=1)
	regressions = regressions.dropna(subset=['pval','neg_log_p']).drop(columns=['note','pval','Conf-interval beta'])

	if model_type == "log":
		regressions.rename(columns={'beta':'LOR'}, inplace=True)
	else:
		regressions.rename(columns={'beta':'Beta'}, inplace=True)

	# add stat info so Explorer interface can disply plots correctly
	regressions.loc[-1, "PheCode"] = "STAT_INFO"
	tmp = dependent + " ~ " + predictors.replace("+", " + ")
	regressions.loc[-1, "Phenotype"] = tmp.replace("PheCode", f"PheCode[{reg_type}]")
	regressions.loc[-1, "Pheno_id"] = model_type

	return regressions





output_columns = [
	'result_type',
	'PheCode',
	'Phenotype',
	'Pheno_id',
	'note',
	'count',
	'neg_log_p',
	'pval',
	'beta',
	'Conf-interval beta',
	'beta_ci_low',
	'beta_ci_up',
	'std_error',
	'category_id',
	'category',
	]

cat_order =  {
	'circulatory system': 0,
	'congenital anomalies': 1,
	'dermatologic': 2,
	'digestive': 3,
	'endocrine/metabolic': 4,
	'genitourinary': 5,
	'hematopoietic': 6,
	'infectious diseases': 7,
	'injuries & poisonings': 8,
	'mental disorders': 9,
	'musculoskeletal': 10,
	'neoplasms': 11,
	'neurological': 12,
	'other': 13,
	'pregnancy complications': 14,
	'respiratory': 15,
	'sense organs': 16,
	'symptoms': 17,
	}

cat_order_df = pd.DataFrame.from_dict(cat_order, orient='index', columns=['category']).reset_index()
cat_order_df.rename(columns={'index':'category_string'}, inplace=True)
phewas_codes.drop(columns='category',inplace=True) # drop old category ID
phewas_codes = phewas_codes.merge(cat_order_df, how='left', on='category_string')






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
