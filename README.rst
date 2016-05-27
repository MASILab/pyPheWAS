pyPhewas
========

pyPhewas module is used to execute a variety of different analyses on large sets of patient phenotype data.

Getting started:

    from pyPhewas import Phewas
   	# define phenotype file, genotype file, or any other desired options
   	p = Phewas(phenotypes, genotypes)
   	p.run_lin() # for a linear regression

Features
--------

- Run linear/logarithmic regressions on patient data
- Plot results

Installation
------------

Install $project by running:

    pip install pyPhewas

Contribute
----------

- Issue Tracker: github.com/BennettLandman/pyPhewas/issues
- Source Code: github.com/BennettLandman/pyPhewas

Support
-------

If you are having issues, please let us know.

License
-------

The project is licensed under the MIT license.
