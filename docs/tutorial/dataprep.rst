Data Preparation
================
This page describes the tools available for preparing your data before running
the pyPheWAS analysis. These tools all require **phenotype** and/or **group**
files. The formats of these files are explained in the :ref:`basics` section.



createGenotypeFile (Create a group file)
----------------------------------------
Split subjects into case (genotype=1) / control (genotype=0) groups based on ICD codes.

Required Arguments:
 * ``--phenotype``: Name of input phenotype file
 * ``--groupout``: Name of output group file
 * ``--case_codes``: Case ICD codes (filename or comma-separated list)
 * ``--code_freq``: Minimum frequency of codes (If 2 comma-separated values are
   given and ctrl_codes is given, 2nd argument is applied to controls)

Optional Arguments [default value]:
 * ``--path``: Path to all input files and destination of output files [current directory]
 * ``--group``: Name of existing group file to add genotype map to
 * ``--ctrl_codes``: Control ICD codes (filename or comma-separated list)

Output:
 * File (*groupout*) containing subject IDs and genotype assignments
 * *Optional:* if an existing file is specified with ``--group``, the genotype
   assignments will be added as a new column to the data in the existing group file.

Specify a list of ICD-9/10 codes that define the case group (genotype=1) and the minimum
frequency of those codes required to be included in the group (e.g. if the
frequency is set to 2, a subject would need to have at least 2 instances of the
case codes in their record to be in the case group). All subjects not in the
case group are put in the control group.

**Example** Define case group as subjects with at least 3 instances of the codes
008 or 134.1; make all other subjects controls::

        createGenotypeFile --phenotype="icd_data.csv" -—groupout="group.csv" --case_codes="008,134.1" --code_freq="3"


Optionally, a list of codes may also be provided for the control group
(genotype=0). In this case, the control group will be composed of subjects not
in the case group that have at least the minimum frequency of control group codes
in their record; *all subjects not in the case or control groups are removed.*
Also optionally, a second argument may be provided to the code frequency input;
if this is specified along with ctrl_codes, the second frequency value will be
applied to the control group.

**Example** Define case group as subjects with at least 3 instances of the codes 008;
define control group as subjects with at least 2 instances of the codes 480.1 or 041::

        createGenotypeFile --phenotype="icd_data.csv" -—groupout="group.csv" --case_codes="008" --ctrl_codes="480.1,041" --code_freq="3,2"


ICD code lists may alternatively be specified by text or csv files. Contents of the
text/csv file should be a comma-separated list similar to the previous examples.
For example, the first example could also be achieved via the following text file and
command:

**case_icd.txt**::

    008,134.1

**Command**::

    createGenotypeFile --phenotype="icd_data.csv" -—groupout="group.csv" --case_codes="case_icd.txt" --code_freq="3"



convertEventToAge (Convert event dates to ages)
-----------------------------------------------
Converts the date at an ICD or CPT event to the subject's age.

Required Arguments:
 * ``--phenotype``:     Name of input phenotype file
 * ``--group``:	        Name of input group file
 * ``--phenotypeout``:  Output phenotype file name (event date replaced with event age)
 * ``--eventcolumn``:	Name of the event date column in the phenotype file
 * ``--etype``:         Type of data (CPT or ICD)

Optional Arguments [default value]:
 * ``--path``:	        The path to all input files and destination of output files [current directory]
 * ``--precision``:	    Decimal precision of the age needed [5]
 * ``--dob_column``:    Name of the date of birth column in the group file ['DOB']

Output:
 * File (phenotypeout) containing ICD/CPT events with event dates replaced by event age

This function converts event dates to event ages, naming the age column according
to the official pyPheWAS phenotype file format ('AgeAtICD' or 'AgeAtCPT').

**Example** Convert CPT event dates to ages with 7-decimal precision::

        convertEventToAge --phenotype="cpt_dates.csv" -—group="group.csv" --phenotypeout="cpt_ages.csv" --eventcolumn="CPT_DATE" --etype="CPT" --precision=7



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
