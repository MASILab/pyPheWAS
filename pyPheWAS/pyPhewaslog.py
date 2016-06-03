# Shikha Chaganti
# Kunal Nabar
# Vanderbilt University
# Medical-image Analysis and Statistical Interpretation Lab

# newphewas
# v2.0

"""
Python implemented code for running logarithmic regressions on PheWAS inputs.

"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import time,math, scipy.stats
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from matplotlib import rcParams


"""
I/O Reading Input From Files
"""
def get_codes():
	"""
	Gets the PheWAS codes from a local csv file and load it into a pandas dataframe.

	"""
	filename = '../resources/codes.csv'	
	return pd.read_csv(filename)

def get_input(path, filename):
	"""
	Read all of the phenotype data from an origin file and load it into a pandas dataframe.
	
	"""
	wholefname = path + filename
	icdfile = pd.read_csv(wholefname)
	g=icdfile.groupby(['id','icd9'])
	idx=g.filter(lambda x: len(x)==1).index
	#icdfile.ix[idx,'icd9']='zzz'
	phenotypes = pd.merge(icdfile,codes,on='icd9')
	return phenotypes

def get_phewas_info(p_index):
	p_code = phewas_codes.loc[p_index].phewas_code	
	corresponding = codes[codes.phewas_code == p_code]

	p_name = corresponding.iloc[0].phewas_string
	p_rollup = ','.join(codes[codes.phewas_code == p_code].icd9.tolist())
	return [p_code, p_name, p_rollup]

def get_group_file(path, filename):
	"""
	Read all of the genotype data from an origin file and load it into a pandas dataframe.
	
	"""
	wholefname = path + filename
	genotypes = pd.read_csv(wholefname)
	return genotypes

"""
Generates a feature matrix of (# of patients)x(icd9 counts)
"""
def generate_feature_matrix(genotypes,phenotypes):
	"""
	Generate the feature matrix from the genotype and phenotype matrix.

	The feature matrix is used to calculate the regressions. 
	The feature matrix has a vector for each ID, each vector has the
	0/1 count as for whether or not rhat phewas code showed up for the given
	patient ID.

	"""
	feature_matrix = np.zeros((genotypes.shape[0],phewas_codes.shape[0]), dtype=int)
	count=0;
	for i in genotypes['id']:
		match=phewas_codes['phewas_code'].isin(list( phenotypes[phenotypes['id']==i]['phewas_code']))
		feature_matrix[count,match[match==True].index]=1
		count+=1
	return feature_matrix

def get_bon_thresh(normalized,power):
	"""
	Calculate the bonferroni correction.

	Divide the power by the sum of all finite values.

	"""
	return power/sum(np.isfinite(normalized))

		

def run_phewas(fm, genotypes,covariates):
	"""
	Calculate an array of p-values.

	Each p-value to the regression performed on each phewas code.
	
	"""
	m = len(fm[0,])
	p_values = np.zeros(m, dtype=float)
	neglogp = np.vectorize(lambda x: -math.log10(x) if x != 0 else 0)
	print('running phewas')
	icodes=[]
	regressions = pd.DataFrame(columns=output_columns)
	for index in range(m):

		phen_vector = fm[:,index]
		
		res=calculate_odds_ratio(genotypes, phen_vector,covariates)
		
		print(index)
		
		phewas_info = get_phewas_info(index)
		stat_info = res[2]
		info = phewas_info[0:2] + stat_info + [phewas_info[2]]
		regressions.loc[index] = info
		
		p_values[index] = res[1]
	normalized = neglogp(p_values)	
	return (np.array(range(m)), normalized, regressions)


"""
Plotting
"""
def get_x_label_positions(categories):
	"""
	Get the location to place each label on the x-axis.

	"""
	tt = Counter(categories)
	s = 0
	label_positions = []
	for _,v in tt.items():
		label_positions.append(s + v//2)
		s += v
	return label_positions

def plot_data_points(x, y, thresh, save=''):
	"""
	Plots x,y, and the threshold line.

	An option can also be made to save the plot in a file of your choosing.

	"""
	c = codes.loc[phewas_codes['index']]
	c = c.reset_index()
	idx = c.sort_values(by='category').index
	x_label_positions = get_x_label_positions(c['category'].tolist())
	x_labels = c.sort_values('category').category_string.drop_duplicates().tolist()
	e = 1
	artists = []
	for i in idx:
		plt.plot(e,y[i],'o', color = plot_colors[c[i:i+1].category_string.values[0]],markersize=10, fillstyle='full', markeredgewidth=0.0)
		if y[i] > thresh:
			artists.append(plt.text(e,y[i],c['phewas_string'][i], rotation=40, va='bottom'))
		e += 1
	plt.axhline(y=-math.log10(0.05), color='blue')
	plt.axhline(y=thresh, color='red')
	plt.xticks(x_label_positions, x_labels,rotation=70, fontsize=10)
	plt.ylim(ymin=0)
	plt.xlim(xmin=0, xmax=len(c))
	plt.ylabel('-log10(p)')
	if save:
		plt.savefig(save,bbox_extra_artists=artists, bbox_inches='tight')
	else:
		plt.subplots_adjust(left=0.05,right=0.85)
		plt.show()
	plt.clf()

def calculate_odds_ratio(genotypes, phen_vector,covariates):
	"""
	Calculates the odds and p-value for using a logarithmic regression.

	"""
	data = genotypes
	data['y']=phen_vector
	f='y~'+covariates
	try:
		logreg = smf.glm(f,data=data,family=sm.families.Binomial()).fit()
		p=logreg.pvalues.genotype
		odds=logreg.deviance	
		conf = logreg.conf_int()
		od = [-math.log10(p), logreg.params.genotype, '[%s,%s]' % (conf[0]['genotype'],conf[1]['genotype'])]
	except:
		odds=0
		p='nan'
		od = ['nan','nan','nan']
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

def phewas(path, filename, groupfile, covariates, save='',output=''):
	"""
	This method will execute a linear regression for phewas.
	
	"""
	start_time = time.time()
	global codes,phewas_codes
	#path,filename,groupfile=('','phenotype.csv','Probable_AD.csv')
	print("reading in data")
	phenotypes = get_input(path, filename)
	genotypes = get_group_file(path, groupfile)
	fm = generate_feature_matrix(genotypes,phenotypes)
	print(len(fm))
	results = run_phewas(fm, genotypes,covariates)
	if output:
		results[2].to_csv(output, index=False)
	thresh = get_bon_thresh(results[1],0.05)
	plot_data_points(results[0],results[1],-math.log10(thresh), save)
	return (results[0], results[1], -math.log10(thresh))


