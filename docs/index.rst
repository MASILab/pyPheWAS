.. pyPheWAS documentation master file, created by
   sphinx-quickstart on Mon May 23 15:05:09 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyPheWAS's documentation!
====================================

The pyPhewas module is used to execute a variety of different analyses on large sets of patient phenotype data.

Getting started::
	
	from pyPhewas import Phewas # imports the Phewas object
	# define phenotype file, genotype file, or any other desired options
	p = Phewas(phenotypes, genotypes)
	p.run_lin() # generates a linear regression for the given data

Features
--------

* Run linear/logarithmic regrssions on patient data
* Plot results

Installation
------------

Install pyPhewas through pip by running::

	pip install pyPhewas

Online
------

* Issue Tracker: github.com/BennettLandman/pyPheWAS/issues
* Source Code: github.com/BennettLandman/pyPheWAS

Support
-------

If you are having issues, please let us know!

License
-------

This project is licensed under the MIT license.

Contents
--------

.. toctree::
   :maxdepth: 2

   tutorial/tutorial
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
* :ref:`modindex`