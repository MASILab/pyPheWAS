Welcome to pyPheWAS's documentation!
====================================

The pyPheWAS module executes PheWAS analyses on large EMR datasets via command
line tools. If this tools contributes to a scientific publication, please cite us.

.. todo:: Add citation.

Features
--------
* Analysis of International Classification of Disease codes (both ICD-9 and ICD-10) and Current Procedural Terminology codes
* EMR data cleaning and preparation
* Compute mass logistic regressions on patient data
* Visualize results
* Examine relative novelty of disease-PheCode associations
* **New!** pyPheWAS Explorer: Interactive visualization of PheDAS experiments
* **New!** Sample synthetic `EMR dataset <https://github.com/MASILab/pyPheWAS/tree/master/synthetic_data>`_

Latest Release: pyPheWAS 4.0.3
------------------------------

This release includes:

* **Novelty Analysis**: examine the relative literary novelty of disease-phecode pairings
* **pyPheWAS Explorer**: an interactive visualization of PheDAS experiments
* createGenotypeFile updated - now called createPhenotypeFile
* Minor bug fixes

Support
-------

* `Issue Tracker <https://github.com/MASILab/pyPheWAS/issues>`_
* `Source Code <https://github.com/MASILab/pyPheWAS>`_

If you are having issues, please let us know! Email me at:

cailey.i.kerley@vanderbilt.edu

License
-------

This project is licensed under the `MIT license <https://github.com/MASILab/pyPheWAS/blob/master/LICENSE>`_.

Contents
--------

.. toctree::
   :maxdepth: 1

   basic
   dataprep
   phewas_tools
   prowas_tools
   novelty_tools
   customizations
   explorer_overview
   references
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
* :ref:`modindex`
