
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
import time,math, scipy.stats
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from matplotlib import rcParams
import os

"""
I/O Reading Input From Files
"""
def get_codes(): #same
	"""
	Gets the PheWAS codes from a local csv file and load it into a pandas DataFrame.

	:returns: All of the codes from the resource file.
	:rtype: pandas DataFrame

	"""
	sep = os.sep
	path = os.path.dirname(os.path.abspath(__file__))
	filename = os.sep.join([path,'..','resources','codes.csv'])
	return pd.read_csv(filename)

def get_input(path, filename): #diff -done - add duration
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

	if  gen_ftype==0:
		g=icdfile.groupby(['id','icd9'])
		idx=g.filter(lambda x: len(x)==1).index
		phenotypes = pd.merge(icdfile,codes,on='icd9')
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
	return phenotypes

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

def generate_feature_matrix(genotypes,phenotypes,control_age): #diff - done
	"""
	Generates the feature matrix that will be used to run the regressions.

	:param genotypes:
	:param phenotypes:
	:type genotypes:
	:type phenotypes:

	:returns:
	:rtype:

	"""
	feature_matrix = np.zeros((genotypes.shape[0],phewas_codes.shape[0]), dtype=int)
	count=0
	for i in genotypes['id']:
		if gen_ftype==0:
			match=phewas_codes['phewas_code'].isin(list( phenotypes[phenotypes['id']==i]['phewas_code']))
			feature_matrix[count,match[match==True].index]=1

		else:
			temp=pd.DataFrame(phenotypes[phenotypes['id']==i][['phewas_code','count','duration']]).drop_duplicates()
			if gen_ftype==1:
				cts = pd.merge(phewas_codes,temp,on='phewas_code',how='left')['count']
				cts[np.isnan(cts)]=0
				feature_matrix[count,:]=cts
			elif gen_ftype==2:
				dura = pd.merge(phewas_codes,temp,on='phewas_code',how='left')['duration']
				dura[np.isnan(dura)]=0
				feature_matrix[count,:]=dura

		count+=1
	return feature_matrix

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


def run_phewas(fm, genotypes ,covariates): #same
	"""
	For each phewas code in the feature matrix, run the specified type of regression and save all of the resulting p-values.

	:param fm: The phewas feature matrix.
	:param genotypes: A pandas DataFrame of the genotype file.
	:param covariates: The covariates that the function is to be run on.

	:returns: A tuple containing indices, p-values, and all the regression data.
	"""
	m = len(fm[0,])
	p_values = np.zeros(m, dtype=float)
	print('running phewas')
	icodes=[]
	# store all of the pertinent data from the regressions
	regressions = pd.DataFrame(columns=output_columns)
	for index in range(m):

		phen_vector = fm[:,index]

		res=calculate_odds_ratio(genotypes, phen_vector,covariates)

		# save all of the regression data
		phewas_info = get_phewas_info(index)
		stat_info = res[2]
		info = phewas_info[0:2] + stat_info + [phewas_info[2]]
		regressions.loc[index] = info

		p_values[index] = res[1]
	return (np.array(range(m)), p_values, regressions)




"""
Plotting
"""
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

def plot_data_points(x, y, thresh0,thresh1,thresh_type, save='', imbalances=np.array([])): #same
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
	plt.axhline(y=-math.log10(0.05), color='black')
	plt.axhline(y=thresh0, color='red')
	plt.axhline(y=thresh1, color='blue')
	plt.xticks(x_label_positions, x_labels, rotation=70, fontsize=10)
	y_label_positions = [-math.log10(0.05), thresh0, thresh1]
	plt.yticks(y_label_positions, ['p=0.05', 'Bonferroni Threshold','Benjamini-Hochberg threshold'], rotation=10,fontsize=10)
	plt.xlim(xmin=0, xmax=len(c))
	plt.ylabel('-log10(p)')

	if thresh_type==0:
		thresh=thresh0
	else:
		thresh=thresh1
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

def calculate_odds_ratio(genotypes, phen_vector,covariates): #diff - done
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
	data['y']=phen_vector
	f='y~'+covariates
	try:
		if gen_ftype==0:
			logreg = smf.glm(f,data=data,family=sm.families.Binomial()).fit()
			p=logreg.pvalues.genotype
			odds=logreg.deviance
			conf = logreg.conf_int()
			od = [-math.log10(p), logreg.params.genotype, '[%s,%s]' % (conf[0]['genotype'],conf[1]['genotype'])]
		else:
			linreg = smf.ols(f,data=data).fit()
			p=linreg.pvalues.genotype
			odds=0
			conf = linreg.conf_int()
			od = [-math.log10(p), linreg.params.genotype, '[%s,%s]' % (conf[0]['genotype'],conf[1]['genotype'])]

	except:
		odds=0
		p=np.nan
		od = [np.nan,np.nan,np.nan]
	return (odds,p,od)

"""
Begin init code
"""

codes = get_codes()
phewas_codes =  pd.DataFrame(codes['phewas_code'].drop_duplicates());
phewas_codes = phewas_codes.reset_index()

output_columns = ['PheWAS Code', 
 'PheWAS Name', 
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

gen_ftype = 0
neglogp = np.vectorize(lambda x: -math.log10(x) if x != 0 else 0)


def phewas(path, filename, groupfile, covariates, reg_type=0, thresh_type=0, control_age=0, save='',output='', show_imbalance=False): #same
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
	:type path: str
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
	global codes,phewas_codes, gen_ftype, neglogp

	print("reading in data")

	gen_ftype = reg_type
	phenotypes = get_input(path, filename)
	genotypes = get_group_file(path, groupfile)
	fm = generate_feature_matrix(genotypes,phenotypes,control_age)
	print(len(fm))
	results = run_phewas(fm, genotypes,covariates)
	regressions = results[2]
	if output:
		regressions.to_csv(output, index=False)
	normalized = neglogp(results[1])
	#if thresh_type==0:
	thresh0 = get_bon_thresh(normalized,0.05)
	#elif thresh_type==1:
	thresh1 = get_fdr_thresh(results[1],0.05)
	imbalances = np.array([])
	if show_imbalance:
		imbalances = get_imbalances(regressions)
	plot_data_points(results[0],normalized,-math.log10(thresh0),-math.log10(thresh1),thresh_type, save, imbalances)
	return (results[0], results[1], imbalances)
