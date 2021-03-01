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

Latest Release: pyPheWAS 3.1.1
------------------------------

This release includes:

* **New Analysis Type:** :ref:`ProWAS Tools`
* **New Plot Type:** Volcano Plot (see :ref:`pyPhewasPlot`)
* :ref:`maximizeControls` now saves explicit Case/Control matches
* New PheCode category colors in plots are more distinguishable
* Improved command line tool argument handling
* Improved error handling
* Documentation overhaul
* API update
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
   tutorial
   dataprep
   phewas_tools
   prowas_tools
   novelty_tools
   references
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
* :ref:`modindex`
