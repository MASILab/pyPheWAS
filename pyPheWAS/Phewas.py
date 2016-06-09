
from pyPhewasv2 import phewas
import pandas as pd
import numpy as np
from pyPhewasv2 import *
import copy

class Phewas:
	
	"""
	This is the object that is used to run Phewas.

	Checks the validity of arguments such as:
	* Valid Phewas file structure
	* Valid values for the arguments passed

	Allows for option such as:
	* Age-matching and control
	* Matching and control by other criteria
	"""
	
	__default_args = {
		'path':'./',
		'covariates':'genotype', 
		'reg_type':0, 
		'thresh_type':0, 
		'save':'',
		'output':'',
		'show_imbalance':False,
	}

	@staticmethod
	def get_default_args():
		return copy.deepcopy(Phewas.__default_args)



	def __verify(self, path, *args):
		path = os.path.abspath(path)
		for element in args:
			floc = os.sep.join([path, element])
			assert os.path.isfile(floc), 'No file found at: %s' % (floc)

	def display_kwargs(self):
		kw = self.kwargs
		pf = self.pfile
		gf = self.gfile

		width = max(len(pf), len(gf), len(kw), 20)
		print('phen'.ljust(width,'.') + pf.rjust(width,'.'))
		print('gen'.ljust(width,'.') + gf.rjust(width,'.'))
		for key in kw:
			left = key.ljust(width,'.')
			right = str(kw[key]).rjust(width, '.')
			print(left + right)

	def set_kwargs(self, **kwargs):
		for key in kwargs:
			self.change_kwarg(key, kwargs[key])

	def change_kwarg(self, key, value):
		kw = self.kwargs
		da = self.__default_args
		assert key in da, '%s is not a valid argument' % (key)
		assert type(kw[key]) == type(da[key]), '%s is not a valid type for %s, use %s' % (type(kw[key]), key, type(da[key]))
		kw[key] = value

	def __init__(self, pfile, gfile, **kwargs):
		"""
		
		
		"""
		self.kwargs = self.get_default_args()
		self.set_kwargs(**kwargs)

		self.__verify(kwargs['path'], pfile, gfile)
		self.pfile = pfile
		self.gfile = gfile

	def run_regression(self):
		self.display_kwargs()
		phewas(self.pfile, self.gfile, **self.kwargs)