#
#
#
#
"""
Phewas object used to run different phenotype statistical anaylses.

This module is used to execute a variety of different analyses on large sets of patient phenotype data.

"""

import pyPhewasv2, pyPhewaslin

class Phewas:
	def __init__(self, inputfile, groupfile, path='', covariates='genotype',save='', output=''):
		self.pargs = [path, inputfile, groupfile, covariates, save, output]
		self.results = None
	def run_lin(self):
		self.result = (pyPhewaslin.phewas(*self.pargs))
	def run_log(self):
		self.result = (pyPhewasv2.phewas(*self.pargs))
	def replot(self):
		if self.result == None:
			print("Please run Phewas before viewing results")
		else:
			pyPhewasv2.plot_data_points(*self.result)