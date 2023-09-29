Welcome to pyPheWAS's documentation!
====================================

The pyPheWAS module executes PheWAS analyses on large EMR datasets via command
line tools. If this tools contributes to a scientific publication, please cite us:

Kerley, C.I., Chaganti, S., Nguyen, T.Q. et al. pyPheWAS: A Phenome-Disease Association Tool for Electronic Medical Record Analysis. *Neuroinformatics* (2022). https://doi.org/10.1007/s12021-021-09553-4

Kerley, C.I., Nguyen T.Q., Ramadass, K, et al. pyPheWAS Explorer: a visualization tool for exploratory analysis of phenome-disease associations. *JAMIA Open* (2023). https://doi.org/10.1093/jamiaopen/ooad018

Features
--------
* Analysis of International Classification of Disease codes (both ICD-9 and ICD-10) and Current Procedural Terminology codes
* EMR data cleaning and preparation
* Compute mass logistic regressions on patient data
* Visualize results
* Examine relative novelty of disease-PheCode associations
* pyPheWAS Explorer: Interactive visualization of PheDAS experiments
* Sample synthetic `EMR dataset <https://github.com/MASILab/pyPheWAS/tree/master/synthetic_data>`_

Latest Release: pyPheWAS 4.2.0
------------------------------

This release includes:

- Default regression equation modified to allow for both canonical and reversed PheWAS equations
- Updated plot styling to improve legibility
- Bug fix: can now run pyPhewasModel/pyProwasModel without covariates
- Other minor bug fixes

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
