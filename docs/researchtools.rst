PheWAS Tools
============
This page describes the command line tools available for running pyPheWAS analyses.
These tools require **phenotype** and **group**
files. The formats of these files are explained in the :ref:`Basics` section.

Using pyPheWAS Research Tools
-----------------------------

The pyPheWAS Research Tools have 3 primary phases. However, age matching and censoring tools exist to change the data beforehand depending on your desired output. The other phases are the Lookup, Model, and Plot phases.

* pyPhewasLookup: generates a feature matrix from the ICD, group data, and type of regression.
* pyPhewasModel: runs the regressions from the feature matrix and outputs a large set of statistical data.
* pyPhewasPlot: this step plots the data from the regressions file.

This is a basic outline of the pyPheWAS research tools structure:

.. figure:: pyPheWAS_Research_Tools.png


Regression Types
----------------
Three regression types are available for the pyPheWAS package.

 1. **log**: binary aggregates (Is a PheCode present for a subject?)
 2. **lin**: count aggregates (How many times is a PheCode present for a subject?)
 3. **dur**: duration aggregates (What is the time interval [years] between the first and last instances of a PheCode present for a subject?)

pyPhewasLookup
--------------

Generate a subject x PheCode feature matrix from ICD data.

Required Arguments:
 * ``--phenotype``: 	Name of the phenotype file
 * ``--group``:		    Name of the group file
 * ``--reg_type``:      Type of regression to use ("log", "lin", or "dur")

Optional Arguments [default value]:
 * ``--path``:		    Path to all input files and destination of output files [current directory]
 * ``--outfile``:	    Base name of the output feature matrix files ["feature_matrix_[group file name].csv"]
 * ``--phewas_cov``:    A PheCodes to use as covariate

Output:
 * Feature matrix with PheCodes as columns and subjects as rows, split into 3 files:

    * **agg_measures**: aggregate PheCode measurement (log/lin/dur)
    * **icd_age**: maximum age on record for each PheCode, may be used as a covariate in pyPhewasModel by specifying "MaxAgeAtICD" in covariates list
    * **phewas_cov**: covarying PheCode matrix, tracks if a subject has at least one record of the PheCode specified by ``phewas_cov``

.. note::
    All 3 feature matrix files are required by pyPhewasModel even if icd_age and
    phewas_cov will be unused.


**Example** Generate a feature matrix for a duration regression::

		pyPhewasLookup  --phenotype="icd_data.csv" --group="group.csv" --reg_type="dur" --outfile="feature_matrix.csv" --path="/Users/me/Documents/EMRdata/"



pyPhewasModel
-------------

pyPhewasModel takes the feature matrices, group file, covariate information, and regression type and runs logarithmic regressions on each PheWAS code to test for association.

Required Arguments:
 * ``--feature_matrix``:the name of the feature matrix CSV file (e.g. "feature_matrix_group.csv")
 * ``--group``:			the name of the group CSV file (e.g. " group.csv")
 * ``--reg_type``:		the regression type to be used ("log", "lin", or "dur")
Optional Arguments:
 * ``--path``:			the path to all input files and destination of output files (*default*: current directory)
 * ``--outfile``:		the name of the output CSV file that contains the regression data (*default*: "regressions_[group file name]")
 * ``--covariates``:	the variables to be used as covariates seperated by '+' (e.g. "SEX" or "SEX+MaxAgeAtICD")
 * ``--response``:	    the variable to predict (instead of genotype)
 * ``--phewas_cov``:	a PheWAS code to use as covariate


A sample execution of *pyPhewasModel*::

		pyPhewasModel --path="/Users/me/Documents/EMRdata/" --feature_matrix="feature_matrix_group.csv" --group="group.csv" --covariates="MaxAgeAtICD" --reg_type="log" --outfile="regressions_group.csv"


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
