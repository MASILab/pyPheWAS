Research Tools
==============

Using pyPheWAS Research Tools
-----------------------------

The pyPheWAS Research Tools have 3 primary phases. However, age matching and censoring tools exist to change the data beforehand depending on your desired output. The other phases are the Lookup, Model, and Plot phases.

* pyPhewasLookup: generates a feature matrix from the icd9, group data, and type of regression.
* pyPhewasModel: runs the regressions from the feature matrix and outputs a large set of statistical data.
* pyPhewasPlot: this step plots the data from the regressions file.

This is a basic outline of the pyPheWAS research tools structure:

.. figure:: pyPheWAS_Research_Tools.png


Regression Types
----------------
Three regression types are available for the pyPheWAS package.

 1. **log**: binary aggregates (Is a PheWAS code present for a subject?)
 2. **lin**: count aggregates (How many times is a PheWAS code present for a subject?)
 3. **dur**: duration aggregates (What is the time interval [years] between the first and last instances of a PheWAS code present for a subject?)

pyPhewasLookup
--------------
 
pyPhewasLookup takes the phenotype file and the group file and will generate the feature matrices.

Required Arguments:
 * ``--phenotype``: 	the name of the phenotype CSV file (e.g. "icd9_data.csv")
 * ``--group``:		    the name of the group CSV file (e.g. "groups.csv")
 * ``--reg_type``:      the type of regression that you would like to use ("log", "lin", or "dur")
Optional Arguments:
 * ``--path``:		    the path to all input files and destination of output files (*default*: current directory)
 * ``--outfile``:	    the name of the output file that contains the feature matrix (*default*: "feature_matrix_[group file name]")
 * ``--phewas_cov``:    a PheWAS code to use as covariate

A sample execution of *pyPhewasLookup*::

		pyPhewasLookup --path="/Users/me/Documents/EMRdata/" --phenotype="icd9_data.csv" --group="group.csv" --reg_type="log" --outfile="feature_matrix_group.csv"


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


