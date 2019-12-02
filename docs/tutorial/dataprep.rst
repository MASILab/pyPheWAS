Data Preparation
================
This page describes the tools available for preparing your data before running the pyPheWAS analysis.

generateGroups (Grouping Tool)
------------------------------

The grouping tool allows you to take two or more icd9 files, and two or more group files. And merge them together, while removing any double counted groups, so that the resulting data files are ready to be run through the pyPheWAS Research Tools.

The options:
 * ``--path``:			        the path to all input files and destination of output files
 * ``--phenotypefiles``:		a list of phenotype file names, each separated by a *+*
 * ``--groupfile``:				a list of group file names, each separated by a *+*
 * ``--phenotypeout``:			the output file name for the merged phenotype files
 * ``--groupout``:				the output file name for the merged group files

A sample execution of *generateGroups*::

		generateGroups --path="/Users/me/Documents/EMRdata" --phenotypefiles="icd9_one.csv+icd9_two.csv" --groupfiles="group_one.csv+group_two.csv" --phenotypeout="new_icd9.csv" --groupout="new_group.csv"

convertEventToAge (Convert event dates to ages)
-----------------------------------------------
Converts event date of ICD9 or CPT to age at the event. Phenotype and group files should be provided with “id” column in both files, and a “DOB” column in the group file.

The options:
 * ``--phenotype``:     phenotype file name
 * ``--group``:	        group file name
 * ``--path``:	        the path to all input files and destination of output files
 * ``--phenotypeout``:  the output file name for the merged phenotype files
 * ``--eventcolumn``:	name of the event date column
 * ``--precision``:	    decimal precision of the age needed
 * ``--type``:          type of data (CPT or ICD)

censorData (ICD/CPT date censoring)
-----------------------------------

Censor files to restrict data to a specific time interval. The default field option is to censor based on AgeAtICD. Can change the default field to other events such as AgeAtDx. 

The options:
 * ``--path``:			the path to all input files and destination of output files
 * ``--phenotype``:		phenotype file name
 * ``--group``:			group file name
 * ``--field``:			the field is the type of event to censor on
 * ``--phenotypeout``:	the output file name for the censored phenotype files
 * ``--groupout``:		the output file name for the censored genotype files
 * ``--start``:			start time for censoring (in years)
 * ``--end``:			end time for censoring (in years)

A sample execution of *censorData*::

		censorData --path="/Users/me/Documents/EMRdata" --phenotype="icd9_data.csv" --group="group.csv" —field=“AgeAtDx” —-phenotypeout="icd9_data_cen.csv" —groupout="group_cen.csv" -—start="0" —-end="2"

maximizeControls (Case/Control matching)
----------------------------------------
Match the subjects in case and control groups based on a criteria such as age ('key'), and on an interval condition ('delta'). The default option for matching groups is genotype (condition='genotype'). The default matching group can be changed to other options such as sex or race.

The options:
 * ``--path``: the path to all input files and destination of output
 * ``--input``:	input group file name
 * ``--output``:	output group file name
 * ``--deltas``:	the intervals for the matching criteria
 * ``--keys``: the fields on which the matching criteria is applied
 * ``--condition``: the field which denotes the groups to be matched
 * ``--goal``: n, indicating the ratio of control and case groups that are being matched
	
A sample execution of * maximizeControls*::

		maximizeControls --path="/Users/me/Documents/EMRdata" --input="group.csv" --output="group__am.csv" --deltas="1,0" --keys="MaxAgeAtVisit+SEX" --condition="genotype" --goal="2"

.. note:: Case/Control matching is performed using the Hopcroft-Karp algorithm. If there are not enough case/control matches, **some case subjects may be dropped**, and will not appear in the output files.

createGenotypeFile (Create a group file)
----------------------------------------
Create a group file by defining ICD-9 codes in the case group and the minimum frequency required to be included in the study.
The options:

 * ``--path``: the path to all input files and destination of output
 * ``--phenotype``: phenotype file name
 * ``--groupout``: output group file name
 * ``--code``: list of ICD-9 codes separated by comma
 * ``--code_freq``: minimum frequency of codes

