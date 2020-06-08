Data Preparation
================
This page describes the command line tools available for preparing your data before running
the pyPheWAS analysis. These tools all require **phenotype** and/or **group**
files. The formats of these files are explained in the :ref:`Basics` section.



createGenotypeFile
------------------
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
 * File (``groupout``) containing subject IDs and genotype assignments
 * *Optional:* if an existing file is specified with ``group``, the genotype
   assignments will be added as a new column to the data in the existing group file.

Specify a list of ICD-9/10 codes that define the case group (genotype=1) and the minimum
frequency of those codes required to be included in the group (e.g. if the
frequency is set to 2, a subject would need to have at least 2 instances of the
case codes in their record to be in the case group). All subjects not in the
case group are put in the control group.

**Example** Define case group as subjects with at least 3 instances of the codes
008 or 134.1; make all other subjects controls::

        createGenotypeFile --case_codes="008,134.1" --code_freq="3" --phenotype="icd_data.csv" -—groupout="group.csv" 


Optionally, a list of codes may also be provided for the control group
(genotype=0) via ``ctrl_codes``. In this case, the control group will be composed of subjects not
in the case group that have at least the minimum frequency of control group codes
in their record; *all subjects not in the case or control groups are removed.*
Also optionally, a second argument may be provided to the ``code_freq`` input;
if this is specified along with ctrl_codes, the second frequency value will be
applied to the control group.

**Example** Define case group as subjects with at least 3 instances of the codes 008;
define control group as subjects with at least 2 instances of the codes 480.1 or 041::

        createGenotypeFile --case_codes="008" --ctrl_codes="480.1,041" --code_freq="3,2" --phenotype="icd_data.csv" -—groupout="group.csv"


ICD code lists may alternatively be specified by text or csv files. Contents of the
text/csv file should be a comma-separated list similar to the previous examples.
For example, the first example could also be achieved via the following text file and
command:

**case_icd.txt**::

    008,134.1

**Command**::

    createGenotypeFile --case_codes="case_icd.txt" --code_freq="3" --phenotype="icd_data.csv" -—groupout="group.csv"



convertEventToAge
-----------------
Converts the date at an ICD or CPT event to the subject's age.

Required Arguments:
 * ``--phenotype``:     Name of input phenotype file
 * ``--group``:	        Name of input group file
 * ``--phenotypeout``:  Name of output phenotype file
 * ``--eventcolumn``:	Name of the event date column in the phenotype file
 * ``--etype``:         Type of data (CPT or ICD)

Optional Arguments [default value]:
 * ``--path``:	        Path to all input files and destination of output files [current directory]
 * ``--precision``:	    Decimal precision of the age needed [5]
 * ``--dob_column``:    Name of the date of birth column in the group file ['DOB']

Output:
 * File (``phenotypeout``) containing ICD/CPT events with event dates replaced by event age

This function converts event dates to event ages, naming the age column according
to the official pyPheWAS phenotype file format ('AgeAtICD' or 'AgeAtCPT').

**Example** Convert CPT event dates to ages with 7-decimal precision::

        convertEventToAge --eventcolumn="CPT_DATE" --etype="CPT" --precision=7 --phenotype="cpt_dates.csv" -—group="group.csv" --phenotypeout="cpt_ages.csv"



censorData
----------
Restrict ICD/CPT data to a specific time interval.

Required Arguments:
 * ``--phenotype``:		Name of phenotype file
 * ``--group``:			Name group file
 * ``--phenotypeout``:	Name of output phenotype file
 * ``--groupout``:		Name of output group file

Optional Arguments [default value]:
 * ``--path``:	        Path to all input files and destination of output files [current directory]
 * ``--efield``:		Name of field in the phenotype file to be censored [AgeAtICD]
 * ``--start``:			Start time for censoring in years [-infinity]
 * ``--end``:			End time for censoring in years [+infinity]
 * ``--delta_field``:	If specified, censor with respect to the interval between delta_field and efield

.. note:: Either the start and/or end arguments must be given.

Output:
 * Phenotype File (``phenotypeout``) containing events censored to specified interval.
 * Group File (``groupout``) containing all subjects with data remaining after censoring.


Specify a range of ages for censoring ICD/CPT event data, such that ``efield`` ages are
censored to the range

        :math:`start \leq efield \leq end`

**Example:** Censor ICD event ages (AgeAtICD) to ages 5 to 18 years-old::

		censorData --start=5 --end=18 --phenotype="icd_data.csv" --group="group.csv" —-phenotypeout="icd_censored.csv" —groupout="group_censored.csv"


Instead of censoring based on absolute age, the user may censor with respect to
another data field using the ``delta_field`` option. If specified, the data is
censored based on the *interval between* ``delta_field`` and ``efield``:

        :math:`start \leq deltafield - efield \leq end`.

**Example:** Censor CPT events to everything previous to 1 year before patient surgery (AgeAtSurgery)::

		censorData --efield="AgeAtCPT" --delta_field="AgeAtSurgery" -—start=1 --phenotype="cpt_data.csv" --group="group.csv" —-phenotypeout="cpt_censored.csv" —groupout="group_censored.csv"



maximizeControls
----------------
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


generateGroups
--------------

The grouping tool allows you to take two or more icd9 files, and two or more group files. And merge them together, while removing any double counted groups, so that the resulting data files are ready to be run through the pyPheWAS Research Tools.

The options:
 * ``--path``:			        the path to all input files and destination of output files
 * ``--phenotypefiles``:		a list of phenotype file names, each separated by a *+*
 * ``--groupfile``:				a list of group file names, each separated by a *+*
 * ``--phenotypeout``:			the output file name for the merged phenotype files
 * ``--groupout``:				the output file name for the merged group files

A sample execution of *generateGroups*::

		generateGroups --path="/Users/me/Documents/EMRdata" --phenotypefiles="icd9_one.csv+icd9_two.csv" --groupfiles="group_one.csv+group_two.csv" --phenotypeout="new_icd9.csv" --groupout="new_group.csv"
