ProWAS Tools
============
This page describes the command line tools available for running ProWAS analyses.
These tools require **phenotype** and **group** files, which are described in the
:ref:`File Formats` section.

Overview of pyProWAS Research Tools
-----------------------------------
Procedural phenome-Wide Association Studies (ProWAS) analyze the relationship between a
large number of EMR phenotypes derived from *Current Procedural Terminology* (CPT)
codes and a dependent variable, typically a genetic marker. It is identical to
PheWAS analysis is all respects, except for 1) the application to CPT codes and 2)
the mapping to procedural phenotype codes (ProCodes). This novel CPT-ProCode
mapping was developed by the MASI lab.

ProWAS consist of four primary phases: 1) data preparation, 2) ProCode mapping
and aggregation, 3) mass ProCode regression, and 4) result visualization. (For more
information, see :ref:`What is PheWAS?`) This page covers phases 2-4, which are
accomplished by the following functions (in order):

* :ref:`pyProwasLookup`: map CPT codes to ProCodes & aggregate
  according to the desired regression type
* :ref:`pyProwasModel`: estimate logistic regression model between dependent variable and
  each ProCode
* :ref:`pyProwasPlot`: visualize the regression results from pyProwasModel

The streamlined tool :ref:`pyProwasPipeline` encompasses all three phases/tools above.

.. note:: For information on the data preparation phase, please see the :ref:`Data Preparation` section.

----------

pyProwasLookup
--------------
Generate a subject x ProCode feature matrix from CPT data.

Maps CPT codes from the phenotype file to their corresponding ProCodes,
then aggregates ProCode data across each subject according to the chosen ``reg_type``.
(Regression types are described in :ref:`Phenotype Aggregation`.)
This is saved as an NxP feature matrix, where N = number of subjects and
P = number of ProCodes.

Required Arguments:
 * ``--phenotype``: 	Name of the phenotype file
 * ``--group``:		    Name of the group file
 * ``--reg_type``:      Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:		    Path to all input files and destination of output files [current directory]
 * ``--outfile``:	    Base name of the output feature matrix files ["feature_matrix _\ ``group``"]
 * ``--prowas_cov``:    A ProCode to use as covariate

Output:
 Feature matrix with ProCodes as columns and subjects as rows, split into 2-3 files

 * **agg_measures**: aggregate ProCode measurement (log/lin/dur)
 * **cpt_age**: maximum age on record for each ProCode, may be used as a covariate
   in **pyProwasModel** by specifying "MaxAgeAtCPT" in covariates list
 * **prowas_cov**: covarying ProCode matrix, tracks if a subject has at least
   one record of the ProCode specified by ``prowas_cov`` (This file will only be
   created if ``prowas_cov`` is provided)


**Example** Generate a feature matrix for a duration regression::

		pyProwasLookup  --reg_type dur --phenotype cpt_data.csv --group group.csv --outfile fm_dur.csv --path /Users/me/Documents/EMRdata/

**Example** Generate a feature matrix for a linear regression with ProCode 67.8 (Laparoscopy) in a covarying feature matrix::

		pyProwasLookup  --reg_type lin --prowas_cov 67.8 --phenotype cpt_data.csv --group group.csv --outfile fm_lin.csv --path /Users/me/Documents/EMRdata/


.. note:: The ``outfile`` argument provides a base name for saving the feature matrix files.
          The three feature matrices are actually saved as
          agg_measures\_\ ``outfile``\ , cpt_age\_\ ``outfile``\ ,
          and prowas_cov\_\ ``outfile``\ .

----------

pyProwasModel
-------------

Perform a mass logistic regression

Iterates over all ProCodes in the feature matrix produced by **pyProwasLookup**
and estimates a regression of the form:

  :math:`procode\_aggregate \sim target + covariates`

or the *reverse* form (`canonical=False`):

  :math:`target \sim procode\_aggregate + covariates`

Linear regression is used if ``reg_type=[lin, dur]`` and ``canonical=True``; otherwise, a logistic regression is used.

.. note:: In version 4.2.0 we changed the default regression equation to the canonical form shown above.
  However, the original pyPheWAS regression equation may still be used via `canonical=False`.

By default, the target variable is 'genotype'; if an alternate variable is specified
by the ``target`` argument, the variable must be a column in the group file.

To use the **cpt_age** feature matrix as a covariate, include 'MaxAgeAtCPT' in
the covariate list. To use the **prowas_cov** feature matrix as a covariate,
specify the ``prowas_cov`` argument. With the exception of these two feature
matrices, all covariates must be columns in the group file.

The saved regression data for each ProCode includes the p-value, -log\ :sub:`10`\ (p-value), beta,
beta's confidence interval, and beta's standard error for the *ProCode_aggregate*
term in the regression model. Additionally, lists of the CPT
codes that map to each ProCode are included.

Regressions are estimated using the [Statsmodels]_ package.

Required Arguments:
 * ``--feature_matrix``: Base name of the feature matrix files
 * ``--group``:			Name of the group file
 * ``--reg_type``:		Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:			Path to all input files and destination of output files [current directory]
 * ``--outfile``:		Name of the output regression data file ["regressions _\ ``group``"]
 * ``--target``:	    Binary variable that indicates case/control groups (default: genotype)
 * ``--covariates``:	Variables to be used as covariates separated by '+' (e.g. "SEX" or "BMI+MaxAgeAtCPT")
 * ``--canonical``:  Use target as a predictor [True, default] or the dependent variable [False] in the regression equation
 * ``--reg_thresh``: Threshold of subjects presenting a ProCode required for running regression (default: 5)
 * ``--prowas_cov``:	A ProCode to use as covariate

Output:
 Regression results for each ProCode saved to the provided ``outfile``

**Example** Compute a duration regression with sex as a covariate::

		pyProwasModel --reg_type dur --covariates sex --feature_matrix fm_dur.csv --group group.csv --outfile regressions_dur.csv --path /Users/me/Documents/EMRdata/

**Example** Compute a binary regression with Dx as the target and sex + cpt_age feature matrix as covariates::

		pyProwasModel --reg_type log --target Dx --covariates sex+MaxAgeAtCPT --feature_matrix my_fm_log.csv --group my_group.csv --outfile reg_log.csv

**Example** Compute a linear regression using the reverse regression equation with the prowas_cov feature matrix for ProCode 67.8 (Laparoscopy) as a covariate::

		pyProwasModel --reg_type lin --prowas_cov 67.8 --canonical False --feature_matrix fm_lin.csv --group my_group.csv --outfile reg_lin_pro678.csv


-----------

pyProwasPlot
------------

Visualizes the regression results through 3 complementary views:

 1. *Manhattan Plot*: This view compares statistical significance across ProCodes.
    ProCodes are presented across the horizontal axis, with -log\ :sub:`10`\ (p) along
    the vertical axis. If ``imbalances = True``\ , marker shape indicates whether
    the effect of each ProCode is positive (+) or negative (-).
 2. *Effect Size Plot*: This view compares effect size across ProCodes. The regression coefficient 
    (or log odds for logistic regressions)
    of each ProCode and its confidence interval are plotted on the horizontal axis,
    with ProCodes presented along the vertical axis. If ``prowas_label = "plot"``\ ,
    ProCode labels are displayed directly on the plot next to their markers. If ``prowas_label = "axis"``\ ,
    ProCodes are displayed outside of the axes, along the left edge.
 3. *Volcano Plot*: This view compares statistical significance and effect size
    across all ProCodes. The effect size of each ProCode is plotted along the
    horizontal axis, with -log\ :sub:`10`\ (p) along the vertical axis.
    ProCodes are colored according to significance level (Not significant, FDR, Bonferroni).

In both the Manhattan and Effect Size plots only ProCodes which are significant
after the chosen multiple comparisons correction is applied are included.

All plots are created using [Matplotlib]_.

Required Arguments:
 * ``--statfile``:		Name of the output regressions file from **pyProwasModel**
 * ``--thresh_type``:	Type of multiple comparisons correction threshold ("bon", "fdr", "custom")

Optional Arguments [default value]:
 * ``--path``:          Path to all input files and destination of output files [current directory]
 * ``--outfile``:       Base name of output plot files [don't save; show interactive plot]
 * ``--imbalance``:     Show the direction of imbalance on the Manhattan plot ([True] or False)
 * ``--plot_all_pts``:  Show all points regardless of significance in the Manhattan plot [True (default) or False]
 * ``--prowas_label``:  Location of the ProCode labels on the Effect Size plot (["plot"] or "axis")
 * ``--old_style``:     Use old plot style (no gridlines, all spines shown)
 * ``--custom_thresh``: Custom threshold value, required if ``thresh_type = "custom"`` (float between 0 and 1)

Threshold Types:
 * *bon*:	    Use the Bonferroni correction
 * *fdr*:	    Use the False Discovery Rate
 * *custom*:	Use a custom threshold specified by ``custom_thresh``

**Example** Plot regression results from the current directory with Bonferroni correction (display results interactively)::

		pyProwasPlot --thresh_type bon --statfile regressions.csv

**Example** Plot regression results with FDR correction and the Log Odds labels displayed on the y-axis (save results)::

		pyProwasPlot --thresh_type fdr --prowas_label axis --outfile my_FDR_plot.eps --statfile regressions.csv --path /Users/me/Documents/EMRdata/

**Example** Plot regression results with a custom threshold and no imbalance on the Manhattan plot (save results)::

		pyProwasPlot --thresh_type custom --custom_thresh 0.001 --imbalance False --outfile my_custom_plot.png --statfile regressions.csv --path /Users/me/Documents/EMRdata/


.. note:: **If outfile is not specified, the plots will not be saved automatically**.
    Instead, all plots will be displayed on the screen by the matplotlib module. It
    is possible to save the plot with any desired file name directly from this display.

.. note:: **Output Formats** Accepted output formats partially depend on which backend is
    active on the user's machine. However, most backends support png, pdf, ps, eps, and svg.
    Vector-based formats (such as svg or svgz) may be opened with image editing software such as Inkscape or
    Photoshop if the user would like to adjust ProCode text locations.

----------

pyProwasPipeline
----------------

**pyProwasPipeline** is a streamlined combination of **pyProwasLookup**, **pyProwasModel**,
and **pyProwasPlot**. If using all default values for optional arguments,
it takes a group file, phenotype file, and regression type and (1) creates the feature
matrix, (2) runs the regressions, and (3) saves Manhattan, Effect Size, and Volcano plots with
both Bonferroni and False Discovery Rate thresholds. All intermediate files
are saved with the ``postfix`` argument appended to the file name.


Required Arguments:
 * ``--phenotype``: 	Name of the phenotype file
 * ``--group``:		    Name of the group file
 * ``--reg_type``:      Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:		    Path to all input files and destination of output files [current directory]
 * ``--postfix``:       Descriptive postfix for output files ["_\ ``covariates``\ _\ ``group``"]
 * ``--target``:	    Binary variable that indicates case/control groups (default: genotype)
 * ``--covariates``:	Variables to be used as covariates separated by '+' (e.g. "SEX" or "BMI+MaxAgeAtCPT")
 * ``--prowas_cov``:    A ProCode to use as covariate
 * ``--canonical``: Use target as a predictor [True, default] or the dependent variable [False] in the regression equation
 * ``--reg_thresh``: Threshold of subjects presenting a ProCode required for running regression (default: 5)
 * ``--thresh_type``:	Type of multiple comparisons correction threshold ("bon", "fdr", "custom")
 * ``--custom_thresh``: Custom threshold value, required if ``thresh_type = "custom"`` (float between 0 and 1)
 * ``--imbalance``:		Show the direction of imbalance on the Manhattan plot ([True] or False)
 * ``--plot_all_pts``: Show all points regardless of significance in the Manhattan plot [True (default) or False]
 * ``--prowas_label``:  Location of the ProCode labels on the Effect Size plot (["plot"] or "axis")
 * ``--old_style``: Use old plot style (no gridlines, all spines shown)
 * ``--plot_format``:   Format for plot files ["png"]


**Example** Run a duration experiment with all default arguments::

		pyProwasPipeline --reg_type dur --phenotype cpt_data.csv --group group.csv

**Example** Run a binary experiment with covariates sex and race, plotting the results with FDR correction, and saving all files with the postfix "binary_prelim"::

		pyProwasPipeline --reg_type log --covariates sex+race --thresh_type fdr --postfix binary_prelim --phenotype cpt_data.csv --group group.csv
