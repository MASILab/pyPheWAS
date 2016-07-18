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

Getting Started
---------------

Install the pyPheWAS package by running::

		pip install pyPheWAS

If this command fails, make sure that you are running Python 2.7+ and that you have pip set up on your machine.

As long as the install is successful, the pyPheWAS package can now be run from any directory.

pyPhewasLookup
--------------
 
pyPhewasLookup takes the phenotype file and the group file and will generate the feature matrix.

The options:
 * ``--path``:		the path to all input files and destination of output files
 * ``--phenotype``:	the name of the phenotype file
 * ``--group``:		the name of the group file
 * ``--outfile``:	the name of the output file that contains the feature matrix
 * ``--reg_type``: the type of regression that you would like to use. (See the key below for more information)

.. note:: the outfile parameter is not required. If it is left off, the default output file will be "feature_matrix_[group file name]" So if "group.csv" is entered as the group file and the outfile parameter is not specified, the feature matrix will be placed in "feature_matrix_group.csv"

The valid options for reg_type:
 * log: logarithmic regression
 * lin: linear regression
 * dur: linear regression on the duration of diseases

A sample execution of *pyPhewasLookup*::

		pyPhewasLookup --path="/Users/me/Documents/EMRdata/" --phenotype="icd9_data.csv" --group="group.csv" --reg_type="log" --outfile="feature_matrix_group.csv"

The "EMRdata" folder before the command:

.. figure:: pyPhewasLookupBefore.png

After the command:

.. figure:: pyPhewasLookupAfter.png

pyPhewasModel
-------------

pyPhewasModel takes the feature matrix, group file, covariate information, and regression type and runs regressions on each PheWAS code to test for association.

The options:
 * ``--path``:			the path to all input files and destination of output files
 * ``--feature_matrix``:the name of the feature matrix file (e.g. "feature_matrix_group.csv")
 * ``--group``:			the name of the group file (e.g. " group.csv")
 * ``--covariates``:	the variables to be used as covariates
 * ``--reg_type``:		the regression type to be used (e.g. "log")
 * ``--outfile``:		the name of the output file that contains the regression data
The valid regression types are listed above under *pyPhewasLookup*

.. note:: the outfile parameter is not required. If it is left off, the default output file will be "regressions_[group file name]" So if "group.csv" is entered as the group file and the outfile parameter is not specified, the feature matrix will be placed in "regressions_group.csv"

.. note:: If multiple covariates are to be used, it is necessary to specify them in one string with a *+* in between them. For example, if you would like to use "genotype" and "age" as covariates, the argument would be ``--covariates="genotype+age"``

A sample execution of *pyPhewasModel*::

		pyPhewasModel --path="/Users/me/Documents/EMRdata/" --feature_matrix="feature_matrix_group.csv" --group="group.csv" --covariates="genotype" --reg_type="log" --outfile="regressions_group.csv"

The "EMRdata" folder before the command:

.. figure:: pyPhewasModelBefore.png

After the command:

.. figure:: pyPhewasLookupAfter.png

pyPhewasPlot
------------

pyPhewasPlot takes the statistics file and threshold type and generates a plot based on the regression data.

The options:
 * ``--path``:			the path to all input files and destination of output files
 * ``--statfile``:		the name of the statistics file
 * ``--imbalance``:		the option of whether or not to show the direction of imbalance in the plot
 * ``--thresh_type``:	the type of threshold to be used in the plot (See the key below for more information)
 * ``--outfile``:		the name of the output file for the plot

.. note:: the outfile is not required. If it is left off, an output file will not be saved to the target directory. Instead, a plot will be displayed on the screen by the matplotlib module. It is possible to save the plot with any desired file name in this display.

.. note:: the ``--imbalance`` option must be either *True* or *False*

The valid options for thresh_type:
 * *bon*:	Use the Bonferroni correction threshold
 * *fdr*:	Use the False Discovery Rate threshold

A sample execution of *pyPhewasPlot*::

		pyPhewasPlot --path="/Users/me/Documents/EMRdata/" --statfile="regressions_group.csv" --imbalance="False" --thresh_type="bon" --outfile="pyPheWAS_plot.png"

The "EMRdata" folder before the command:

.. figure:: pyPhewasPlotBefore.png

After the command:

.. figure:: pyPhewasPlotAfter.png


Additional Research Tools
=========================

Grouping Tool (generateGroups)
-------------

The grouping tool allows you to take two or more icd9 files, and two or more group files. And merge them together, while removing any double counted groups, so that the resulting data files are ready to be run through the pyPheWAS Research Tools.

The options:
 * ``--path``:			the path to all input files and destination of output files
 * ``--phenotypefiles``:		a list of phenotype file names, each separated by a *+*
 * ``--groupfiles``:				a list of group file names, each separated by a *+*
 * ``--phenotypeout``:			the output file name for the merged phenotype files
 * ``--groupout``:				the output file name for the merged group files

A sample execution of *generateGroups*::

		generateGroups --path="/Users/me/Documents/EMRdata" --phenotypefiles="icd9_one.csv+icd9_two.csv" --groupfiles="group_one.csv+group_two.csv" --phenotypeout="new_icd9.csv" --groupout="new_group.csv"

Age Matching
------------



Censoring
---------

The age censoring tool allows you to censor your data by age. The tool allows you to do two things:

#. Censor the data from exact start and end ages
#. Censor the data as a range from a period of time, relative to an age in the demographic data.

..note:: The field 'AgeAtICD' must be included in the phenotype file, this is the age that will be censored according to the parameters.

The options:
 * ``--path``:			the path to all input files and destination of the output files
 * ``--phenotype``:		the name of the phenotype file
 * ``--group``:			the name of the group file
 * ``--phenotypeout``:	the output file name for the set of age censored data points
 * ``--start``:			the start time of the range
 * ``--end``:			the end time of the range
 * ``--field``:			the name of the field from which the difference will be calculated (*optional*, leave blank if using exact ages)

Exact start and end ages:
	This method *does not* require the ``field`` input. ``start`` will the beginning of the range and ``end`` will be the end of the age range. For example if ``--start=2`` and ``--end=4``, then only data that has 'AgeAtICD' between the 2 and 4 will be included in the output field.

A sample execution of *censorData* without using ``field``::

	censorData --path="/Users/me/Documents/EMRdata/" --phenotype="icd9_data.csv" --group="group.csv" --start=2 --end=4 --phenotypeout="newicd9.csv"

Delta from the ``field`` parameter:
	This method *does* require the ``field`` input. The difference between the ``field`` parameter and the 'AgeAtICD' parameter will be calculated for each data point. Then only data points for which the difference falls in between start and end will be included. For example, if the parameter ``field=AgeAtDx``, ``start=2`` and ``end=4``. Only data points for which 'AgeAtICD' is 2 to 4 years prior to 'AgeAtDx' will be included. Suppose for a given patient, 'AgeAtDx'=10, for that patient, only data points from the icd9 file in which 'AgeAtICD' are between 6 and 8 will be included.

A sample execution of *censorData* with using ``field``::

	censorData --path="/Users/me/Documents/EMRdata/" --phenotype="icd9_data.csv" --group="group.csv" --start=2 --end=4 --phenotypeout="newicd9.csv" --field="AgeAtDx"

..note:: In order to use the second censoring method. The parameter ``field`` must be included as a column in the group file, and must be an age, similar to the 'AgeAtICD'

Event to Age
------------

A lot of EMR data comes in the form of event dates instead of ages for each set of data. To amend this, this tool exists to convert all of the event dates to ages (*AgeAtICD*). It takes the difference from DOB and event date and fills the age for each data point.

..note:: Date of birth must be included in the group data, referred to as *DOB*.

The options:
 * ``--path``:			the path to all input files and destination of output files
 * ``--phenotype``:		the name of the phenotype file
 * ``--group``:			the name of the group file
 * ``--phenotypeout``:	the name of the output file
 * ``--eventcolumn``:	the name of the column in the phenotype file that includes event dates
 * ``--precision``:		the number of decimal points to be included in the age

A sample execution of *convertEventToAge*::

		convertEventToAge --path="/Users/me/Documents/EMRData/ --phenotype="icd9_data.csv" --group="group.csv" --phenotypeout="icd9_with_age.csv" --eventcolumn="Event_date" --precision=2

