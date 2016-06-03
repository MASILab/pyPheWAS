#
#
#
#
"""
Phewas object used to run different phenotype statistical anaylses.

This module is used to execute a variety of different functio

"""

import pyPhewasv2, pyPhewaslin

class Phewas:

	"""
	Phewas object, used to execute different types of phenotype analyses. 
	
	"""

	def __init__(self, inputfile, groupfile, path='', covariates='genotype',save='', output=''):
		"""
		Initializes the Phewas object so that different analyses can be run.

		Takes in required Phewas arguments such as an inputfile and groupfile
		as well as a variety of other arguments such as path to files, covariates
		an output destination for the plot created, and an output destination for the 
		regression.

		"""
		self.pargs = [path, inputfile, groupfile, covariates, save, output]
		self.results = None
	def run_lin(self):
		"""
		Runs a linear regression using the current Phewas arguments.

		"""
		self.result = (pyPhewaslin.phewas(*self.pargs))
	def run_log(self):
		"""
		Runs a logarithmic regression using the current Phewas arguments.

		"""
		self.result = (pyPhewasv2.phewas(*self.pargs))
	def replot(self):
		"""
		Replots the most recent analysis that has been run.

		Some of the Phewas analyses can be time consuming. This prevents the user
		from having to re-execute the analysis every time they close a graph and would
		like to view it again.
		
		"""
		if self.result == None:
			print("Please run Phewas before viewing results")
		else:
			pyPhewasv2.plot_data_points(*self.result)