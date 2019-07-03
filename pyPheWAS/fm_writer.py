"""
Feature Matrix Writer Class
====================================
Author: Cailey Kerley
Date: July 2, 2019

Class for writing the pyPheWAS feature matrices
"""
import os
import csv


class FM_Writer:
	def __init__(self, path, outfile,key_list):
		agg_out = os.path.join(path, 'agg_measures_' + outfile)
		self.agg_file = open(agg_out, 'w')
		self.agg_writer = csv.DictWriter(self.agg_file, key_list)
		icd_age_out = os.path.join(path, 'icd_age_' + outfile)
		self.icd_age_file = open(icd_age_out, 'w')
		self.icd_age_writer = csv.DictWriter(self.icd_age_file, key_list)
		phewas_cov_out = os.path.join(path, 'phewas_cov_' + outfile)
		self.phewas_cov_file = open(phewas_cov_out, 'w')
		self.phewas_cov_writer = csv.DictWriter(self.phewas_cov_file, key_list)


	def print_fm(self, agg, icd_age, phewas_cov):
		"""
		Prints feature matrix data for one subject.
		:param agg: aggregate feature matrix
		:param icd_age: age at ICD feature matrix
		:param phewas_cov: co-varying phewas code feature matrix
		:type agg: dictionary
		:type icd_age: dictionary
		:type phewas_cov: dictionary

		:returns: None
		"""

		self.agg_writer.writerows(agg)
		self.icd_age_writer.writerows(icd_age)
		self.phewas_cov_writer.writerows(phewas_cov)

		return

	def __del__(self):
		self.agg_file.close()
		self.icd_age_file.close()
		self.phewas_cov_file.close()