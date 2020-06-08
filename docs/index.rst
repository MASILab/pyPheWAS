.. pyPheWAS documentation master file, created by
   sphinx-quickstart on Mon May 23 15:05:09 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyPheWAS's documentation!
====================================

The pyPhewas module is used to execute a variety of different analyses on large sets of patient phenotype data.

Features
--------
* Data cleaning and preparation
* Run logarithmic regressions on patient data
* Plot results

Latest Release: pyPheWAS 3.0.0
------------------------------

This release includes:

* Support for both ICD 9 and ICD 10
* All 3 regression types (binary, count, & duration) optimized for big data
* pyPhewasPipeline: a streamlined combination of pyPhewasLookup, pyPhewasModel, and pyPhewasPlot
* Compatibility with Python 3
* Age matching now saves the explicit mapping between controls/cases in addition to the resulting group file
* Operation of the ICD censoring function matches the description in the documentation
* Minor bug fixes

Support
-------

* Issue Tracker: github.com/MASILab/pyPheWAS/issues
* Source Code: github.com/MASILab/pyPheWAS

If you are having issues, please let us know! Email me at:

cailey.i.kerley@vanderbilt.edu

License
-------

This project is licensed under the MIT license.

Contents
--------

.. toctree::
   :maxdepth: 1

   basic
   dataprep
   researchtools
   tutorial
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
* :ref:`modindex`
