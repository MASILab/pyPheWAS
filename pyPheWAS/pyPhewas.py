"""
Kunal Nabar
Vanderbilt University
Medical-image Analysis and Statistical Interpretation Lab

newphewas
v0.1

"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import time,math, scipy.stats
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

"""
I/O Reading Input From Files
"""
def get_codes():
	filename = 'codes.csv'
	f = open(filename)	
	data = []
	f.readline()
	for line in f:
		data.append(line.split(',')[0:5])
	return zip(*data)

def get_input(path, filename):
	f = open(path + filename)

	data = []
	f.readline()
	for line in f:
		data.append(line.strip().split(',')[0:3])
	p_id, p_icd9, p_event_count = [list(x) for x in zip(*data)]
	p_event_count = [int(x) for x in p_event_count]

	associations = [p_id, p_icd9, p_event_count]
	return associations

def get_group_file(path, filename):
	f = open(path + filename)

	data = []
	f.readline()
	for line in f:
		data.append(line.strip().split(','))
	group_id, group_bool = [list(x) for x in zip(*data)]
	group_bool = np.array(group_bool).astype(int)
	group_data = [group_id, group_bool]
	return group_id, group_bool

"""

"""

def icd9_to_name(icd9):
	return icd9_string[icd9_index_table[icd9]]

"""
Create a mapping from icd9 to index
"""
def generate_icd9_index_table():
	icd9_index_table = {}
	for index, code in enumerate(unique_icd9):
		icd9_index_table[code] = index
	return icd9_index_table

def generate_patient_index_table():
	patient_index_table = {}
	for index,patient_id in enumerate(unique_id):
		patient_index_table[patient_id] = index
	return patient_index_table
"""
Create a lookup table that maps id to a list of tuples of order
(icd9, count)
"""
def generate_lookup_table(p_id, p_icd9, p_event_count):
	patient_icd9_lookup_table = defaultdict(lambda: [])
	for index, code in enumerate(p_id):
		if p_icd9[index].find('.') == 2:
			p_icd9[index] = '0' + p_icd9[index]
		patient_icd9_lookup_table[code].append((p_icd9[index],p_event_count[index]))
	return patient_icd9_lookup_table

"""
Create the feature vector for a given id
Uses the the icd9 lookup table to assign the values of the feature vector
"""
def generate_feature_vector(identification):
	feature_vector = np.zeros(NUM_UNIQUE_ICD9, dtype=int)
	for code, count in patient_icd9_lookup_table[identification]:
		try:
			index = icd9_index_table[code]
		except KeyError:
			continue
		feature_vector[index] = count
	return feature_vector

"""
Generates a feature matrix of (# of patients)x(icd9 counts)
"""
def generate_feature_matrix():
	feature_matrix = np.zeros((len(unique_id),NUM_UNIQUE_ICD9), dtype=int)
	for index, element in enumerate(unique_id):
		feature_matrix[index,] = generate_feature_vector(element)
	return feature_matrix

def get_boolean_vector(fv):
	return (fv > 0).astype(int)


def run_phewas(fm, group_id, group_bool):
	m = len(fm[0,])
	p_values = np.zeros(m, dtype=float)
	neglogp = np.vectorize(lambda x: -math.log10(x) if x != 0 else 0)
	print('running phewas')
	icodes=[]
	for index in range(m):
		#print(index)
		icd9_vector = fm[:,index]
		#res=calculate_odds_ratio(group_id, group_bool, icd9_vector)
		#p_values[index] = res[1]
		#icodes.append(unique_icd9[index])	
	normalized = neglogp(p_values)
	return (np.array(range(m)), normalized)

"""
Plotting
"""


def plot_data_points(x, y, point_style='ro'):
	print('generating plots')
	print(max(y))
	plt.plot(x, y, point_style)
	for i,code in enumerate(unique_icd9):
		if y[i] > 1.3:
			plt.text(x[i],y[i],icd9_to_name(code),rotation=40,va='bottom')
	plt.show()


"""
odds ratio
			group_bool (f1)
    			0     1
 			0   A  |  B
icd9 (f2) 	  -----|-----
 			1   C  |  D
"""

def nCr(n,r):
	f = math.factorial
	try:
		return f(n)//(f(r)*f(n-r))
	except:
		print(n,r)

def myFisher(a,b,c,d):
	n = a+b+c+d
	numerator = nCr(a+b,a) * nCr(c+d,c)
	denominator = nCr(n,a+c)
	return (numerator/denominator)

def calculate_odds_ratio(group_id, group_bool, icd9_vector):
	# A corresponds to not having the icd9 code or being in the group
	# B corresponds to not having the icd9 code but being in the group
	# C corresponds to having the icd9 but not being in the group
	# D corresponds to having the icd9 and being in the group
	#
	# D/C => having both / having just the icd9
	# B/A => having autism and not the icd9 / no icd9 no autism
	# the 1st corresponds to the ratio of the population who has the icd9 code also has the condition
	# the 2nd corresponds to the ratio of the population who has the condition but not the icd9
	icd9_vector = get_boolean_vector(icd9_vector)
	#print(len(grop_bool),len(icd9_vector))
	#f1 = np.where(group_bool == 1)[0]
	#f2 = np.where(icd9_vector == 1)[0]
	#f3 = np.where(group_bool == 0)[0]
	#f4 = np.where(icd9_vector == 0)[0]
	#A = len(np.intersect1d(f3, f4))
	#B = len(np.intersect1d(f1, f4))
	#C = len(np.intersect1d(f2, f3))
	#D = len(np.intersect1d(f1, f2))
	odds = 0
	data = pd.Series([icd9_vector,group_bool],index=['y','X'])
	data=data.dropna()


	try:
		#odds,p,dof,ex = scipy.stats.chi2_contingency([[A, B], [C, D]])
	#	odds,p = scipy.stats.fisher_exact([[A,B],[C,D]])
		logreg = smf.glm('y~X',data=data,family=sm.families.Binomial()).fit()
		p=logreg.pvalues.X
		odds=logreg.deviance
	except:
		odds=0
		p=1
	if np.isnan(p):
		p=1
	return (odds,p)

"""
Begin init code
"""

unique_icd9, icd9_string, phewas_code, desc, rollup = [list(x) for x in get_codes()]
unique_id = []
patient_index_table = {}
patient_icd9_lookup_table = {}
associations = []
icd9_index_table = generate_icd9_index_table()
NUM_UNIQUE_ICD9 = len(unique_icd9)

def phewas(path, filename, groupfile):
	# the path and filename of the goal file 
	# must hold the following format
	# patient id, icd9 code, event count
	# filename = 'testdata.csv'
	start_time = time.time()
	global unique_id,associations,patient_index_table,patient_icd9_lookup_table
	#path,filename,groupfile=('','phenotype.csv','Probable_AD.csv')
	print("reading in data")
	associations = get_input(path, filename)
	group_id,group_bool = get_group_file(path, groupfile)
	unique_id = sorted(list(set(group_id)))
	print("generating lookup table")
	patient_index_table = generate_patient_index_table()
	patient_icd9_lookup_table = generate_lookup_table(*associations)
	print("generating feature matrix")
	fm = generate_feature_matrix()
	print(len(fm))
	#results = run_phewas(fm, group_id, group_bool)

	#print('--- %s seconds --- ' % (time.time() - start_time))

	#plot_data_points(*results)



