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
def get_codes(): #same
	"""
	Gets the PheWAS codes from a local csv file and load it into a pandas dataframe.
	"""
	filename = '../resources/codes.csv'	
	return pd.read_csv(filename)


def get_input(path, filename): #diff -done - add duration
	"""
	Read all of the phenotype data from an origin file and load it into a pandas dataframe.
	
	"""
	wholefname = path + filename
	icdfile = pd.read_csv(wholefname)
	
	if  gen_ftype==0:
		g=icdfile.groupby(['id','icd9'])
		idx=g.filter(lambda x: len(x)==1).index
		phenotypes = pd.merge(icdfile,codes,on='icd9')
	else:
		phenotypes = pd.merge(icdfile,codes,on='icd9')	
		phenotypes['count']=0
		phenotypes['count']=phenotypes.groupby(['id','phewas_code'])['count'].transform('count')					
		phenotypes['duration']=phenotypes.groupby(['id','phewas_code'])['AgeAtICD'].transform('max')-phenotypes.groupby(['id','phewas_code'])['AgeAtICD'].transform('min')+1
	return phenotypes

def get_phewas_info(p_index): #same
	
	p_code = phewas_codes.loc[p_index].phewas_code	
	corresponding = codes[codes.phewas_code == p_code]

	p_name = corresponding.iloc[0].phewas_string
	p_rollup = ','.join(codes[codes.phewas_code == p_code].icd9.tolist())
	return [p_code, p_name, p_rollup]

def get_group_file(path, filename): #same
	wholefname = path + filename
	genotypes = pd.read_csv(wholefname)
	return genotypes

def get_imbalances(regressions):
	"""
	Generates a numpy array of the imbalances.
	
	The mapping from the beta of the regressions to the imbalance is as follows
	For any given value *x*
	*x* == nan -> 0
	*x* < 0 -> -1
	*x* > 0 -> +1
	
	These are then used in conjunction with imbalance_colors in the plotting method to show imbalance via colors.
	"""

	imbalance = np.array(regressions['beta'])
	imbalance[np.isnan(imbalance)] = 0
	imbalance[imbalance > 0] = 1
	imbalance[imbalance < 0] = -1
	return imbalance

"""
Generates a feature matrix of (# of patients)x(icd9 counts)
"""
def generate_feature_matrix(genotypes,phenotypes): #diff - done
	feature_matrix = np.zeros((genotypes.shape[0],phewas_codes.shape[0]), dtype=int)
	count=0;
	for i in genotypes['id']:
		if gen_ftype==0:
			match=phewas_codes['phewas_code'].isin(list( phenotypes[phenotypes['id']==i]['phewas_code']))
			feature_matrix[count,match[match==True].index]=1
			
		else:
			temp=pd.DataFrame(phenotypes[phenotypes['id']==i][['phewas_code','count']]).drop_duplicates()
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
	return power/sum(np.isfinite(normalized))


def get_fdr_thresh(p_values, power):
	sn = np.sort(p_values)
	sn = sn[np.isfinite(sn)]
	sn = sn[::-1]
	for i in range(len(sn)):
		thresh=0.05*i/len(sn)
		if sn[i]<=power:
			break
	return sn[i]
		

def run_phewas(fm, genotypes,covariates): #same
	m = len(fm[0,])
	p_values = np.zeros(m, dtype=float)
	print('running phewas')
	icodes=[]
	regressions = pd.DataFrame(columns=output_columns)
	for index in range(m):

		phen_vector = fm[:,index]
		
		res=calculate_odds_ratio(genotypes, phen_vector,covariates)
		
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

def plot_data_points(x, y, thresh, save, imbalances=[]): #same
	c = codes.loc[phewas_codes['index']]
	c = c.reset_index()
	idx = c.sort_values(by='category').index
	x_label_positions = get_x_label_positions(c['category'].tolist(), imbalances.size == 0)
	x_labels = c.sort_values('category').category_string.drop_duplicates().tolist()
	e = 1
	artists = []
	for i in idx:
		if imbalances.size == 0:
			plt.plot(e,y[i],'o', color=plot_colors[c[i:i+1].category_string.values[0]],markersize=10, fillstyle='full', markeredgewidth=0.0)
		else:
			plt.plot(e,y[i], 'o', color=imbalance_colors[imbalances[i]], fillstyle='full', markeredgewidth=0.0)
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

def calculate_odds_ratio(genotypes, phen_vector,covariates): #diff - done
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
	1: 'mediumspringgreen',
	-1: 'red'
}

gen_ftype = 0
neglogp = np.vectorize(lambda x: -math.log10(x) if x != 0 else 0)


def phewas(path, filename, groupfile, covariates, reg_type=0, thresh_type=0, save='',output='', show_imbalance=False): #same
	# the path and filename of the goal file 
	# must hold the following format
	# patient id, icd9 code, event count
	# filename = 'testdata.csv'
	start_time = time.time()
	global codes,phewas_codes, gen_ftype, neglogp
	
	print("reading in data")

	gen_ftype = reg_type
	phenotypes = get_input(path, filename)
	genotypes = get_group_file(path, groupfile)
	fm = generate_feature_matrix(genotypes,phenotypes)
	print(len(fm))
	results = run_phewas(fm, genotypes,covariates)
	if output:
		results[2].to_csv(output, index=False)
	normalized = neglogp(results[1])	
	if thresh_type==0:
		thresh = get_bon_thresh(normalized,0.05)
	elif thresh_type==1:
		thresh = get_fdr_thresh(results[1],0.05)
	if show_imbalance:
		imbalances = get_imbalances(results[2])
	plot_data_points(results[0],normalized,-math.log10(thresh), save)
	return (results[0], results[1], -math.log10(thresh))


