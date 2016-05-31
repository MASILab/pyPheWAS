pyPhewas
========

pyPhewas module is used to execute a variety of different analyses on large sets of patient phenotype data.

Getting started::

	from pyPhewas import Phewas # imports the Phewas object
	# define phenotype file, genotype file, or any other desired options
	p = Phewas(phenotypes, genotypes)
	p.run_lin() # generates a linear regression for the given data

Features
--------

* Run linear/logarithmic regressions on patient data
* Plot results

Installation
------------

Install pyPhewas by running::

    pip install pyPhewas

Contribute
----------

* Issue Tracker: github.com/BennettLandman/pyPheWAS/issues
* Source Code: github.com/BennettLandman/pyPheWAS

Support
-------

If you are having issues, please let us know.

License
-------

The project is licensed under the MIT license.
