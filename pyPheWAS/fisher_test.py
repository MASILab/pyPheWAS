import math
import scipy.stats
import numpy as np

def get_cont_table(group_id,group_bool,icd9_vector):
	#icd9_vector = get_boolean_vector(icd9_vector)
	f1 = np.where(group_bool == 1)[0]
	f2 = np.where(icd9_vector == 1)[0]
	f3 = np.where(group_bool == 0)[0]
	f4 = np.where(icd9_vector == 0)[0]
	A = len(np.intersect1d(f3, f4))
	B = len(np.intersect1d(f1, f4))
	C = len(np.intersect1d(f2, f3))
	D = len(np.intersect1d(f1, f2))
	return (A,B,C,D)

def get_data(fm, group_bool, group_id):
	x = []
	y = []
	for i in range(0,len(fm[1])):
		print(i)
		if sum(fm[:,i]) == 0:
			continue
		A,B,C,D = get_cont_table(group_id, group_bool, fm[:,i])
		u,v = getBoth(A,B,C,D)
		if u != 1.0:
			x.append(u)
			y.append(v)
	return(x,y)

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

def theirFisher(a,b,c,d):
	(odds,p) = scipy.stats.fisher_exact([[a,b],[c,d]])
	return p


def getBoth(a,b,c,d):
	return (theirFisher(a,b,c,d),myFisher(a,b,c,d))