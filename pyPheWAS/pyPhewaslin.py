"""
Shikha Chaganti
Kunal Nabar
Vanderbilt University
Medical-image Analysis and Statistical Interpretation Lab

newphewas
v2.0

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
	filename = 'codes.csv'	
	return pd.read_csv(filename)

def get_input(path, filename):
	wholefname = path + filename
	icdfile = pd.read_csv(wholefname)
	g=icdfile.groupby(['id','icd9'])
	phenotypes = pd.merge(icdfile,codes,on='icd9')	
	phenotypes['count']=0
	phenotypes['count']=phenotypes.groupby(['id','phewas_code'])['count'].transform('count')
	return phenotypes

def get_group_file(path, filename):
	wholefname = path + filename
	genotypes = pd.read_csv(wholefname)
	return genotypes

"""
Generates a feature matrix of (# of patients)x(icd9 counts)
"""
def generate_feature_matrix(genotypes,phenotypes):
	feature_matrix = np.zeros((genotypes.shape[0],phewas_codes.shape[0]), dtype=int)
	count=0;
	for i in genotypes['id']:
		temp=pd.DataFrame(phenotypes[phenotypes['id']==i][['phewas_code','count']]).drop_duplicates()
		cts = pd.merge(phewas_codes,temp,on='phewas_code',how='left')['count']
		cts[np.isnan(cts)]=0
		feature_matrix[count,:]=cts
		count+=1
	return feature_matrix

def get_bon_thresh(normalized,power):
	return power/sum(np.isfinite(normalized))



def run_phewas(fm, genotypes,covariates):
	m = len(fm[0,])
	p_values = np.zeros(m, dtype=float)
	neglogp = np.vectorize(lambda x: -math.log10(x) if x != 0 else 0)
	print('running phewas')
	icodes=[]
	for index in range(m):

		phen_vector = fm[:,index]
		
		res=calculate_odds_ratio(genotypes, phen_vector,covariates)
		
		print(index)
		
			
		
		p_values[index] = res[1]
	normalized = neglogp(p_values)	
	return (np.array(range(m)), normalized)

"""
Plotting
"""
def get_x_label_positions(categories):
	tt = Counter(categories)
	s = 0
	label_positions = []
	for _,v in tt.items():
		print(v//2)
		label_positions.append(s + v//2)
		s += v
	return label_positions

def plot_data_points(x, y, thresh, save):
	c = codes.loc[phewas_codes['index']]
	c = c.reset_index()
	idx = c.sort_values(by='category').index
	x_label_positions = get_x_label_positions(c['category'].tolist())
	x_labels = c.sort_values('category').category_string.drop_duplicates().tolist()
	e = 1
	for i in idx:
		plt.plot(e,y[i],'o', color = plot_colors[c[i:i+1].category_string.values[0]],markersize=10, fillstyle='full', markeredgewidth=0.0)
		if y[i] > thresh:
			plt.text(e,y[i],c['phewas_string'][i], rotation=40, va='bottom')
		e += 1
	plt.axhline(y=-math.log10(0.05), color='blue')
	plt.axhline(y=thresh, color='red')
	plt.xticks(x_label_positions, x_labels,rotation=70, fontsize=15)
	plt.ylim(ymin=0)
	plt.xlim(xmin=0, xmax=len(c))
	plt.ylabel('-log10(p)')
	rcParams.update({'figure.autolayout': True})
	#plt.tight_layout(pad=0.2,w_pad=1,h_pad=0)
	if save != None:
		plt.savefig(save)
	plt.show()


def calculate_odds_ratio(genotypes, phen_vector,covariates):
	
	data = genotypes
	data['y']=phen_vector
	f='y~'+covariates
	try:
		linreg = smf.ols(f,data=data).fit()
		p=linreg.pvalues.genotype
		odds=0
	
	except:
		odds=0
		p='nan'
	return (odds,p)

"""
Begin init code
"""

codes = get_codes()
phewas_codes =  pd.DataFrame(codes['phewas_code'].drop_duplicates());
phewas_codes = phewas_codes.reset_index()
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

def phewas(path, filename, groupfile,covariates,save=''):
	# the path and filename of the goal file 
	# must hold the following format
	# patient id, icd9 code, event count
	# filename = 'testdata.csv'
	start_time = time.time()
	global codes,phewas_codes
	#path,filename,groupfile=('','phenotype.csv','Probable_AD.csv')
	print("reading in data")
	phenotypes = get_input(path, filename)
	genotypes = get_group_file(path, groupfile)
	fm = generate_feature_matrix(genotypes,phenotypes)
	print(len(fm))
	results = run_phewas(fm, genotypes,covariates)
	thresh = get_bon_thresh(results[1],0.05)
	plot_data_points(results[0],results[1],-math.log10(thresh),save)



