PheWAS Tools
============
This page describes the command line tools available for running PheWAS analyses.
These tools require **phenotype** and **group** files, which are described in the
:ref:`File Formats` section.

Overview of pyPheWAS Research Tools
-----------------------------------
Phenome-Wide Association Studies (PheWAS) compare EMR phenotypes with a single dependent variable,
historically a genetic marker. The tools described on this page are specifically for the study of
phenotypes derived from *International Classification of Disease* (ICD) codes.
They employ a mapping from ICD-9 and ICD-10 codes to 1,866 hierarchical phenotype codes (PheCodes).
This mapping was originally constructed solely for ICD-9 codes [Denny2013]_,
with later improvements to the ICD-9 mapping [Wei2017]_ and the addition of the
ICD-10 code mapping [Wu2019]_. More information about the mapping can be found
at the `PheWAS Catalog website <https://phewascatalog.org>`_.

PheWAS consist of four primary phases: 1) data preparation, 2) PheCode mapping
and aggregation, 3) mass PheCode regression, and 4) result visualization. (For more
information, see :ref:`What is PheWAS?`) This page
covers phases 2-4, which are accomplished by the following functions (in order):

* :ref:`pyPhewasLookup`: map ICD-9 and ICD-10 codes to PheCodes & aggregate
  according to the desired regression type
* :ref:`pyPhewasModel`: estimate logistic regression model between genotype and
  each PheCode
* :ref:`pyPhewasPlot`: visualize the regression results from pyPhewasModel

The streamlined tool :ref:`pyPhewasPipeline` encompasses all three phases/tools above.

.. note:: For information on the data preparation phase, please see the :ref:`Data Preparation` section.


pyPhewasLookup
--------------
Generate a subject x PheCode feature matrix from ICD data.

Maps ICD-9 and ICD-10 codes from the phenotype file to their corresponding PheCodes,
then aggregates PheCode data across each subject according to the chosen ``reg_type``.
(Regression types are described in :ref:`Phenotype Aggregation`.)
This is saved as an NxP feature matrix, where N = number of subjects and
P = number of PheCodes.

The ICD-9 map currently being used is `version 1.2 <https://phewascatalog.org/phecodes>`_
and the ICD-10 map is `version 1.2beta <https://phewascatalog.org/phecodes_icd10>`_.

Required Arguments:
 * ``--phenotype``: 	Name of the phenotype file
 * ``--group``:		    Name of the group file
 * ``--reg_type``:      Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:		    Path to all input files and destination of output files [current directory]
 * ``--outfile``:	    Base name of the output feature matrix files ["feature_matrix _\ ``group``"]
 * ``--phewas_cov``:    A PheCode to use as covariate

Output:
 Feature matrix with PheCodes as columns and subjects as rows, split into 2-3 files

 * **agg_measures**: aggregate PheCode measurement (log/lin/dur)
 * **icd_age**: maximum age on record for each PheCode, may be used as a covariate
   in **pyPhewasModel** by specifying "MaxAgeAtICD" in covariates list
 * **phewas_cov**: covarying PheCode matrix, tracks if a subject has at least one
   record of the PheCode specified by ``phewas_cov`` (This file will only be
   created if ``phewas_cov`` is provided)


**Example** Generate a feature matrix for a duration regression::

		pyPhewasLookup  --reg_type="dur" --phenotype="icd_data.csv" --group="group.csv" --outfile="fm_dur.csv" --path="/Users/me/Documents/EMRdata/"

**Example** Generate a feature matrix for a linear regression with PheCode 495 (Asthma) in a covarying feature matrix::

		pyPhewasLookup  --reg_type="lin" --phewas_cov="495" --phenotype="icd_data.csv" --group="group.csv" --outfile="fm_lin.csv" --path="/Users/me/Documents/EMRdata/"


.. note:: The ``outfile`` argument provides a base name for saving the feature matrix files.
          The three feature matrices are actually saved as
          agg_measures\_\ ``outfile``\ , icd_age\_\ ``outfile``\ ,
          and phewas_cov\_\ ``outfile``\ .


pyPhewasModel
-------------

Perform a mass logistic regression

Iterates over all PheCodes in the feature matrix produced by **pyPhewasLookup**
and estimates a logistic regression of the form:

    :math:`Pr(response) \sim logit(PheCode\_aggregate + covariates)`

By default, the response variable is 'genotype'; if an alternate variable is specified
by the ``response`` argument, the variable must be a column in the group file.

To use the **icd_age** feature matrix as a covariate, include 'MaxAgeAtICD' in
the covariate list. To use the **phewas_cov** feature matrix as a covariate,
specify the ``phewas_cov`` parameter. With the exception of these two feature
matrices, all covariates must be included as columns in the group file.

The saved regression data for each PheCode includes the p-value, -log\ :sub:`10`\ (p-value), beta,
beta's confidence interval, and beta's standard error for the *PheCode_aggregate*
term in the logit model. Additionally, lists of the ICD-9/ICD-10
codes that map to each PheCode are included.

Logistic regressions are estimated using the [Statsmodels]_ package.

Required Arguments:
 * ``--feature_matrix``: Base name of the feature matrix files
 * ``--group``:			Name of the group file
 * ``--reg_type``:		Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:			Path to all input files and destination of output files [current directory]
 * ``--outfile``:		Name of the output regression data file ["regressions _\ ``group``"]
 * ``--response``:	    Variable to predict ['genotype']
 * ``--covariates``:	Variables to be used as covariates separated by '+' (e.g. "SEX" or "BMI+MaxAgeAtICD")
 * ``--phewas_cov``:	A PheCode to use as covariate

Output:
 Regression results for each PheCode saved to the provided ``outfile``

**Example** Compute a duration regression with sex as a covariate::

		pyPhewasModel --reg_type="dur" --covariates="sex" --feature_matrix="fm_dur.csv" --group="group.csv" --outfile="regressions_dur.csv" --path="/Users/me/Documents/EMRdata/"

**Example** Compute a binary regression with sex and the icd_age feature matrix as covariates::

		pyPhewasModel --reg_type="log" --covariates="sex+MaxAgeAtICD" --feature_matrix="my_fm_log.csv" --group="my_group.csv" --outfile="reg_log.csv"

**Example** Compute a linear regression with the phewas_cov feature matrix for PheCode 495 (Asthma) as a covariate::

		pyPhewasModel --reg_type="lin" --phewas_cov="495" --feature_matrix="fm_lin.csv" --group="my_group.csv" --outfile="reg_lin_phe495.csv"


.. note:: To prevent false positives & improve statistical power, regressions
          are only computed for PheCodes which present in greater than 5
          subjects. PheCodes which do not meet this criteria are
          not included in the output regression file.

.. note:: For phenotypes that present in both the case (``response`` = 1) and
          control (``response`` = 0) groups, maximum likelihood optimization is
          used to compute the logistic regression. For phenotypes that only
          present in one of those groups, regularized maximum likelihood
          optimization is used.


pyPhewasPlot
------------

Visualizes the regression results through 3 complementary views:

1. *Manhattan Plot*: This view compares statistical significance across PheCodes.
   PheCodes are presented across the horizontal axis, with -log\ :sub:`10`\ (p) along
   the vertical axis. If ``imbalances = True``\ , marker shape indicates whether
   the effect of each PheCode is positive (+) or negative (-).
2. *Log Odds Plot*: This view compares effect size across PheCodes. The log odds
   of each PheCode and its confidence interval are plotted on the horizontal axis,
   with PheCodes presented along the vertical axis. If ``phewas_label = "plot"``\ ,
   PheCode labels are displayed directly on the plot next to their markers. If ``phewas_label = "axis"``\ ,
   PheCodes are displayed outside of the axes, along the left edge.
3. *Volcano Plot*: This view compares statistical significance and effect size
   across all PheCodes. The log odds of each PheCode is plotted along the
   horizontal axis, with -log\ :sub:`10`\ (p) along the vertical axis.
   PheCodes are colored according to significance level (Not significant, FDR, Bonferroni).

In both the Manhattan and Log Odds plots:

* PheCode markers are colored and sorted according to 18 general categories
  (mostly organ systems and disease groups, e.g. “circulatory system” and
  “mental disorders”).
* Only PheCodes which are significant after the chosen multiple comparisons
  correction is applied are included.

All plots are created using [Matplotlib]_.

Required Arguments:
 * ``--statfile``:		Name of the output regressions file from **pyPhewasModel**
 * ``--thresh_type``:	Type of multiple comparisons correction threshold ("bon", "fdr", "custom")

Optional Arguments [default value]:
 * ``--path``:          Path to all input files and destination of output files [current directory]
 * ``--outfile``:       Base name of output plot files [don't save; show interactive plot]
 * ``--imbalance``:		Show the direction of imbalance on the Manhattan plot ([True] or False)
 * ``--phewas_label``:  Location of the PheCode labels on the Log Odds plot (["plot"] or "axis")
 * ``--custom_thresh``: Custom threshold value, required if ``thresh_type = "custom"`` (float between 0 and 1)

Threshold Types:
 * *bon*:	    Use the Bonferroni correction
 * *fdr*:	    Use the False Discovery Rate
 * *custom*:	Use a custom threshold specified by ``custom_thresh``

**Example** Plot regression results from the current directory with Bonferroni correction (display results interactively)::

		pyPhewasPlot --thresh_type="bon" --statfile="regressions.csv"

**Example** Plot regression results with FDR correction and the Log Odds labels displayed on the y-axis (save results)::

		pyPhewasPlot --thresh_type="fdr" --phewas_label="axis" --outfile="my_FDR_plot.eps" --statfile="regressions.csv" --path="/Users/me/Documents/EMRdata/"

**Example** Plot regression results with a custom threshold and no imbalance on the Manhattan plot (save results)::

		pyPhewasPlot --thresh_type="custom" --custom_thresh=0.001 --imbalance=False --outfile="my_custom_plot.png" --statfile="regressions.csv" --path="/Users/me/Documents/EMRdata/"


.. note:: **If outfile is not specified, the plots will not be saved automatically**.
    Instead, all plots will be displayed on the screen by the matplotlib module. It
    is possible to save the plot with any desired file name directly from this display.

.. note:: **Output Formats** Accepted output formats partially depend on which backend is
    active on the user's machine. However, most backends support png, pdf, ps, eps, and svg.
    Vector-based formats (such as svg or svgz) may be opened with image editing software such as Inkscape or
    Photoshop if the user would like to adjust PheCode text locations.

pyPhewasPipeline
----------------

**pyPhewasPipeline** is a streamlined combination of **pyPhewasLookup**, **pyPhewasModel**,
and **pyPhewasPlot**. If using all default values for optional arguments,
it takes a group file, phenotype file, and regression type and (1) creates the feature
matrix, (2) runs the regressions, and (3) saves Manhattan, Log Odds, and Volcano plots with
both Bonferroni and False Discovery Rate thresholds. All intermediate files
are saved with the ``postfix`` argument appended to the file name.


Required Arguments:
 * ``--phenotype``: 	Name of the phenotype file
 * ``--group``:		    Name of the group file
 * ``--reg_type``:      Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:		    Path to all input files and destination of output files [current directory]
 * ``--postfix``:       Descriptive postfix for output files ["_\ ``covariates``\ _\ ``group``"]
 * ``--response``:	    Variable to predict ['genotype']
 * ``--covariates``:	Variables to be used as covariates separated by '+' (e.g. "SEX" or "BMI+MaxAgeAtICD")
 * ``--phewas_cov``:    A PheCode to use as covariate
 * ``--thresh_type``:	Type of multiple comparisons correction threshold ("bon", "fdr", "custom")
 * ``--imbalance``:		Show the direction of imbalance on the Manhattan plot ([True] or False)
 * ``--phewas_label``:  Location of the PheCode labels on the Log Odds plot (["plot"] or "axis")
 * ``--custom_thresh``: Custom threshold value, required if ``thresh_type = "custom"`` (float between 0 and 1)
 * ``--plot_format``:   Format for plot files ["png"]


**Example** Run a duration experiment with all default arguments::

		pyPhewasPipeline --reg_type="dur" --phenotype="icd_data.csv" --group="group.csv"

**Example** Run a binary experiment with covariates sex and race, plotting the results with FDR correction, and saving all files with the postfix "binary_prelim"::

		pyPhewasPipeline --reg_type="log" --covariates="sex+race" --thresh_type="fdr" --postfix="binary_prelim" --phenotype="icd_data.csv" --group="group.csv"
