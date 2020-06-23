PheWAS Tools
============
This page describes the command line tools available for running PheWAS analyses.
These tools require **phenotype** and **group**
files. The formats of these files are explained in the :ref:`Basics` section.

Using pyPheWAS Research Tools
-----------------------------

PheWAS analysis consists of 4 primary phases (illustrated below): 1) data preparation, 2) PheCode mapping
and aggregation, 3) mass PheCode regression, and 4) result visualization. This page
covers phases 2-4, which are accomplished by the following functions (in order):

* **pyPhewasLookup**: map ICD-9 and ICD-10 codes to PheCodes & aggregate
  according to the desired regression type
* **pyPhewasModel**: estimate logistic regression model between genotype and
  each PheCode
* **pyPhewasPlot**: visualize the regression results from pyPhewasModel


.. figure:: pyPheWAS_Research_Tools.png


.. note:: For information on the data preparation phase, please see the :ref:`Data Preparation` section.



Regression Type
---------------
Three regression types are available for PheWAS analyses.

 1. **log**: binary aggregates (Is a PheCode present for a subject?)
 2. **lin**: count aggregates (How many times is a PheCode present for a subject?)
 3. **dur**: duration aggregates (What is the time interval [years] between the first and last instances of a PheCode for a subject?)


pyPhewasLookup
--------------
Generate a subject x PheCode feature matrix from ICD data.

Maps ICD-9 and ICD-10 codes from the phenotype file to their corresponding PheCodes,
then aggregates PheCode data across each subject according to the chosen :ref:`Regression Type`.
This is saved as a 3xNxP feature matrix, where N = number of subjects and
P = number of PheCodes.

Required Arguments:
 * ``--phenotype``: 	Name of the phenotype file
 * ``--group``:		    Name of the group file
 * ``--reg_type``:      Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:		    Path to all input files and destination of output files [current directory]
 * ``--outfile``:	    Base name of the output feature matrix files [*feature_matrix_[group file name].csv*]
 * ``--phewas_cov``:    A PheCode to use as covariate

Output:
 Feature matrix with PheCodes as columns and subjects as rows, split into 3 files

 * **agg_measures**: aggregate PheCode measurement (log/lin/dur)
 * **icd_age**: maximum age on record for each PheCode, may be used as a covariate in pyPhewasModel by specifying "MaxAgeAtICD" in covariates list
 * **phewas_cov**: covarying PheCode matrix, tracks if a subject has at least one record of the PheCode specified by ``phewas_cov``


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

Iterates over all PheCodes in the feature matrix produced by ``pyPhewasLookup``
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

Required Arguments:
 * ``--feature_matrix``: Base name of the feature matrix files
 * ``--group``:			Name of the group file
 * ``--reg_type``:		Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:			Path to all input files and destination of output files [current directory]
 * ``--outfile``:		Name of the output regression data file [*regressions_[group file name].csv*]
 * ``--response``:	    Variable to predict ['genotype']
 * ``--covariates``:	Variables to be used as covariates separated by '+' (e.g. "SEX" or "BMI+MaxAgeAtICD")
 * ``--phewas_cov``:	a PheWAS code to use as covariate

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

pyPhewasPlot takes the regressions file and threshold type and generates two plots (Manhattan and Log Odds) based on the regression data.

Required Arguments:
 * ``--statfile``:		the name of the regressions CSV file
 * ``--imbalance``:		whether or not to show the direction of imbalance in the plot ("True" or "False")
 * ``--thresh_type``:	the type of threshold to be used in the plot (See the key below for more information)
Optional Arguments:
 * ``--path``:          the path to all input files and destination of output files (*default*: current directory)
 * ``--outfile``:       the name of the output PDF file for the plot
 * ``--custom_thresh``: custom threshold value, required if *thresh_type* ="custom" (float between 0 and 1)


The valid options for thresh_type:
 * *bon*:	    Use the Bonferroni correction threshold
 * *fdr*:	    Use the False Discovery Rate threshold
 * *custom*:	Use a custom threshold specified by *--custom_thresh*

.. note:: **If outfile is not specified, the plot will not be saved automatically**. Instead, a plot will be displayed on the screen by the matplotlib module. It is possible to save the plot with any desired file name in this display.


A sample execution of *pyPhewasPlot*::

		pyPhewasPlot --path="/Users/me/Documents/EMRdata/" --statfile="regressions_group.csv" --imbalance="False" --thresh_type="bon" --outfile="pyPheWAS_plot.pdf"

pyPhewasPipeline
----------------

pyPhewasPipeline is a streamlined combination of pyPhewasLookup, pyPhewasModel, and pyPhewasPlot. If using all default
values for the optional arguments, it takes a group file, phenotype file, and regression type and (1) creates the feature
matrix, (2) runs the regressions, and (3) saves Manhattan and Log Odds plots with both the BonFerroni and False Discovery
Rate thresholds. All intermediate files are saved with the *postfix* argument appended to the file name.


Required Arguments:
 * ``--phenotype``: 	the name of the phenotype CSV file (e.g. "icd9_data.csv")
 * ``--group``:			the name of the group CSV file (e.g. " group.csv")
 * ``--reg_type``:		the regression type to be used ("log", "lin", or "dur")
Optional Arguments:
 * ``--path``:          the path to all input files and destination of output files (*default*: current directory)
 * ``--postfix``:       descriptive postfix for output files (*default*: "[covariates]_[group file name]")
 * ``--phewas_cov``:    a PheWAS code to use as covariate
 * ``--covariates``:	the variables to be used as covariates seperated by '+' (e.g. "SEX" or "SEX+MaxAgeAtICD")
 * ``--response``:	    the variable to predict (instead of genotype)
 * ``--imbalance``:		whether or not to show the direction of imbalance in the plot, must be "True" or "False" (*default*: True)
 * ``--thresh_type``:	the type of threshold to be used in the plot (See the key below for more information)
 * ``--custom_thresh``: custom threshold value, required if *thresh_type* ="custom" (float between 0 and 1)


The valid options for thresh_type:
 * *bon*:	    Use the Bonferroni correction threshold
 * *fdr*:	    Use the False Discovery Rate threshold
 * *custom*:	Use a custom threshold specified by *--custom_thresh*


A sample execution of *pyPhewasPlot*::

		pyPhewasPipline --path="/Users/me/Documents/EMRdata/" --phenotype="icd9_data.csv" --group="group.csv" --reg_type="log" --postfix="poster_Nov22"
