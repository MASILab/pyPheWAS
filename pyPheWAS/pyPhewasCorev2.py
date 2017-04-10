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
import sys


def get_codes(): #same
	"""
	Gets the PheWAS codes from a local csv file and load it into a pandas DataFrame.

	:returns: All of the codes from the resource file.
	:rtype: pandas DataFrame

	"""
	sep = os.sep
	path = os.path.dirname(os.path.abspath(__file__))
	filename = os.sep.join([path,'resources','codes.csv'])
	return pd.read_csv(filename)

def get_group_file(path, filename): #same
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


def get_input(path, filename, reg_type): #diff -done - add duration
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
	if reg_type == 0:
		phenotypes = pd.merge(icdfile,codes,on='icd9')
		phenotypes['MaxAgeAtICD'] = 0
		phenotypes['MaxAgeAtICD'] = phenotypes.groupby(['id', 'phewas_code'])['AgeAtICD'].transform('max')
	else:
		"""
		This needs to be changed, need to adjust for a variety of different naming conventions
		in the phenotype file, not simply 'AgeAtICD', 'id', 'icd9', etc. 
		Either we need to adjust for different names in the code, or state explicitly in the 
		documentation that we cannot do things like this.
		"""
		phenotypes = pd.merge(icdfile,codes,on='icd9')	
		phenotypes['count']=0
		phenotypes['count']=phenotypes.groupby(['id','phewas_code'])['count'].transform('count')					
		phenotypes['duration']=phenotypes.groupby(['id','phewas_code'])['AgeAtICD'].transform('max')-phenotypes.groupby(['id','phewas_code'])['AgeAtICD'].transform('min')+1
		phenotypes['MaxAgeAtICD'] = 0
		phenotypes['MaxAgeAtICD'] = phenotypes.groupby(['id', 'phewas_code'])['AgeAtICD'].transform('max')
	return phenotypes

def generate_feature_matrix(genotypes,phenotypes, reg_type,phewas_cov=''): #diff - done
	"""
	Generates the feature matrix that will be used to run the regressions.

	:param genotypes:
	:param phenotypes:
	:type genotypes:
	:type phenotypes:

	:returns: 
	:rtype:

	"""
	feature_matrix = np.zeros((3,genotypes.shape[0],phewas_codes.shape[0]), dtype=int)
	count = 0;
	for i in genotypes['id']:
		if reg_type == 0:
			temp=pd.DataFrame(phenotypes[phenotypes['id']==i][['phewas_code','MaxAgeAtICD']]).drop_duplicates()
			match=phewas_codes['phewas_code'].isin(list( phenotypes[phenotypes['id']==i]['phewas_code']))
			feature_matrix[0][count,match[match==True].index]=1
			age = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['MaxAgeAtICD']
			age[np.isnan(age)] = genotypes[genotypes['id'] == i].iloc[0]['MaxAgeAtVisit']
			feature_matrix[1][count, :] = age
			if phewas_cov:
				feature_matrix[2][count, :] = int(phewas_cov in list(phenotypes[phenotypes['id'] == i]['phewas_code']))

		else:
			if reg_type == 1:
				temp=pd.DataFrame(phenotypes[phenotypes['id']==i][['phewas_code','MaxAgeAtICD','count']]).drop_duplicates()
				cts = pd.merge(phewas_codes,temp,on='phewas_code',how='left')['count']
				cts[np.isnan(cts)]=0
				feature_matrix[0][count,:]=cts
				age = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['MaxAgeAtICD']
				age[np.isnan(age)] = genotypes[genotypes['id']==i].iloc[0]['MaxAgeAtVisit']
				feature_matrix[1][count, :] = age
				if phewas_cov:
					feature_matrix[2][count, :] = int(
						phewas_cov in list(phenotypes[phenotypes['id'] == i]['phewas_code']))

			elif reg_type==2:
				temp=pd.DataFrame(phenotypes[phenotypes['id']==i][['phewas_code','MaxAgeAtICD','duration']]).drop_duplicates()
				dura = pd.merge(phewas_codes,temp,on='phewas_code',how='left')['duration']
				dura[np.isnan(dura)]=0
				feature_matrix[0][count,:]=dura
				age = pd.merge(phewas_codes, temp, on='phewas_code', how='left')['MaxAgeAtICD']
				age[np.isnan(age)] = genotypes[genotypes['id']==i].iloc[0]['MaxAgeAtVisit']
				feature_matrix[1][count, :] = age
				if phewas_cov:
					feature_matrix[2][count, :] = int(
						phewas_cov in list(phenotypes[phenotypes['id'] == i]['phewas_code']))

		count+=1
	return feature_matrix

"""

Statistical Modeling

"""
def get_phewas_info(p_index): #same
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

def calculate_odds_ratio(genotypes, phen_vector1,phen_vector2,reg_type,covariates,response='',phen_vector3=''): #diff - done
	"""
	Runs the regression for a specific phenotype vector relative to the genotype data and covariates.

	:param genotypes: a DataFrame containing the genotype information
	:param phen_vector: a array containing the phenotype vector
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
	data['y']=phen_vector1
	data['MaxAgeAtICD'] = phen_vector2
	#f='y~'+covariates
	if response:
		f = response+'~ y + ' + covariates
		if phen_vector3.any():
			data['phe'] = phen_vector3
			f = response + '~ y + phe +' + covariates
	else:
		f = 'y ~' + covariates
		if phen_vector3.any():
			data['phe'] = phen_vector3
			f = 'y ~ phe +' + covariates
	try:
		if reg_type==0:
			logreg = smf.logit(f,data).fit(method='bfgs',disp=False)
			p=logreg.pvalues.genotype
			odds=logreg.params.genotype
			conf = logreg.conf_int()
			od = [-math.log10(p), logreg.params.genotype, '[%s,%s]' % (conf[0]['genotype'],conf[1]['genotype'])]
		else:
			linreg = smf.logit(f,data).fit(method='bfgs',disp=False)
			p=linreg.pvalues.genotype
			odds=0
			conf = linreg.conf_int()
			od = [-math.log10(p), linreg.params.genotype, '[%s,%s]' % (conf[0]['genotype'],conf[1]['genotype'])]
	except:
		odds=0
		p=np.nan
		od = [np.nan,np.nan,np.nan]
	return (odds,p,od)

def run_phewas(fm, genotypes ,covariates, reg_type, response='',phewas_cov=''): #same
	"""
	For each phewas code in the feature matrix, run the specified type of regression and save all of the resulting p-values.
	
	:param fm: The phewas feature matrix.
	:param genotypes: A pandas DataFrame of the genotype file.
	:param covariates: The covariates that the function is to be run on.

	:returns: A tuple containing indices, p-values, and all the regression data.
	"""
	m = len(fm[0,0])
	p_values = np.zeros(m, dtype=float)
	icodes=[]
	# store all of the pertinent data from the regressions
	regressions = pd.DataFrame(columns=output_columns)
	for index in range(m):
		print index
		phen_vector1 = fm[0][:,index]
		phen_vector2 = fm[1][:,index]
		phen_vector3 = fm[2][:, index]
		res=calculate_odds_ratio(genotypes, phen_vector1,phen_vector2, reg_type,covariates,response=response,phen_vector3=phen_vector3)

		# save all of the regression data
		phewas_info = get_phewas_info(index)
		stat_info = res[2]
		info = phewas_info[0:2] + [res[1]] + stat_info  + [phewas_info[2]]
		regressions.loc[index] = info
		
		p_values[index] = res[1]
	return regressions

def get_bon_thresh(normalized,power): #same
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
	return power/sum(np.isfinite(normalized))


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
		thresh=0.05*i/len(sn)
		if sn[i]<=power:
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

def get_x_label_positions(categories, lines=True): #same
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
	for _,v in tt.items():
		if lines:
			inc = v//2
		else:
			inc = v
		label_positions.append(s + inc)
		s += v
	return label_positions

def plot_data_points(y, thresh, save='', imbalances=np.array([])): #same
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
	show_imbalance = imbalances.size != 0
	
	# Sort the phewas codes by category.
	c = codes.loc[phewas_codes['index']]
	c = c.reset_index()
	idx = c.sort_values(by='category').index

	# Get the position of the lines and of the labels
	linepos = get_x_label_positions(c['category'].tolist(), True)
	x_label_positions = get_x_label_positions(c['category'].tolist(), False)
	x_labels = c.sort_values('category').category_string.drop_duplicates().tolist()
	
	# Plot each of the points, if necessary, label the points.
	e = 1
	artists = []
	for i in idx:
		if show_imbalance:
			plt.plot(e,y[i], 'o', color=imbalance_colors[imbalances[i]], fillstyle='full', markeredgewidth=0.0)
		else:
			plt.plot(e,y[i],'o', color=plot_colors[c[i:i+1].category_string.values[0]],markersize=10, fillstyle='full', markeredgewidth=0.0)
		if y[i] > thresh:
			artists.append(plt.text(e,y[i],c['phewas_string'][i], rotation=40, va='bottom'))
		e += 1

	# If the imbalance is to be shown, draw lines to show the categories.
	if show_imbalance:
		for pos in linepos:
			plt.axvline(x=pos, color='black', ls='dotted')

	# Plot a blue line at p=0.05 and plot a red line at the line for the threshold type.
	plt.axhline(y=-math.log10(0.05), color='blue')
	plt.axhline(y=thresh, color='red')

	# Set windows and labels
	plt.xticks(x_label_positions, x_labels,rotation=70, fontsize=10)
	plt.ylim(ymin=0)
	plt.xlim(xmin=0, xmax=len(c))
	plt.ylabel('-log10(p)')

	# Determine the type of output desired (saved to a plot or displayed on the screen)
	if save:
		pdf = PdfPages(save)
		pdf.savefig(bbox_extra_artists=artists, bbox_inches='tight')
		pdf.close()
	else:
		plt.subplots_adjust(left=0.05,right=0.85)
		plt.show()

	# Clear the plot in case another plot is to be made.
	plt.clf()

def process_args(kwargs,optargs,*args):
	clean = np.vectorize(lambda x: x[x.rfind('-')+1:] + '=')
	searchfor = clean(list(optargs.keys()))
	opts, rem = getopt.getopt(args, '',searchfor)
	assert len(rem) == 0, 'Unknown arguments included %s' % (str(rem))
	for option in opts:
		k,v = option
		kwargs[optargs[k]] = v

	return kwargs

def display_kwargs(kwargs):
	print ("Arguments: ")
	for k,v in kwargs.items():
		left = str(k).ljust(30,'.')
		right = str(v).rjust(50,'.')
		print(left + right)

output_columns = ['PheWAS Code', 
 'PheWAS Name',
 'p-val',
 '\"-log(p)\"',
 'beta',
 'Conf-interval beta',
 'ICD-9']
plot_colors = {'-' : 'gold',
 'circulatory system' : 'red',
 'congenital anomalies': 'mediumspringgreen',
 'dermatologic' : 'maroon',
 'digestive' : 'green',
 'endocrine/metabolic' : 'darkred',
 'genitourinary' : 'black',
 'hematopoietic' : 'orange',
 'infectious diseases' : 'blue',
 'injuries & poisonings' : 'slategray',
 'mental disorders' : 'fuchsia',
 'musculoskeletal' : 'darkgreen',
 'neoplasms' : 'teal',
 'neurological' : 'midnightblue',
 'pregnancy complications' : 'gold',
 'respiratory' : 'brown',
 'sense organs' : 'darkviolet',
 'symptoms' : 'darkviolet'}
imbalance_colors = {
	0: 'white',
	1: 'deepskyblue',
	-1: 'red'
}
regression_map = {
	'log':0,
	'lin':1,
	'lind':2
}
threshold_map = {
	'bon': 0,
	'fdr': 1
}

codes = get_codes()
phewas_codes =  pd.DataFrame(codes['phewas_code'].drop_duplicates());
phewas_codes = phewas_codes.reset_index()


