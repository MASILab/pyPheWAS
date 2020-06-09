Data Preparation
================
This page describes the command line tools available for preparing your data before running
a PheWAS or ProWAS analysis. These tools all require **phenotype** and/or **group**
files. The formats of these files are explained in the :ref:`Basics` section.


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

**Example** Censor ICD event ages (AgeAtICD) to ages 5 to 18 years-old::

		censorData --start=5 --end=18 --phenotype="icd_data.csv" --group="group.csv" —-phenotypeout="icd_censored.csv" —groupout="group_censored.csv"


Instead of censoring based on absolute age, the user may censor with respect to
another data field using the ``delta_field`` option. If specified, the data is
censored based on the *interval between* ``delta_field`` and ``efield``:

        :math:`start \leq deltafield - efield \leq end`.

**Example** Censor CPT events to everything previous to 1 year before patient surgery (AgeAtSurgery)::

		censorData --efield="AgeAtCPT" --delta_field="AgeAtSurgery" -—start=1 --phenotype="cpt_data.csv" --group="group.csv" —-phenotypeout="cpt_censored.csv" —groupout="group_censored.csv"


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



maximizeControls
----------------
Match subjects in case and control groups based on group variables.

Required Arguments:
 * ``--input``:     Name of input group file
 * ``--keys``:      Comma-separated list of matching criteria (must be columns in group file)
 * ``--deltas``:	Comma-separated list of tolerance intervals for the matching criteria
 * ``--goal``:      n, target matching ratio (control:case => n:1)

Optional Arguments [default value]:
 * ``--path``:      Path to all input files and destination of output files [current directory]
 * ``--output``:	Name of output group file [input__matched.csv]
 * ``--condition``: Field denoting group assignments [genotype]

Output:
 * Group file (``output``) containing only matched cases/controls.
 * Match file (output__matched_pairs.csv) containing explicit case to control match mapping.

Match cases/controls based on similarity in matching criteria via the Hopcroft-Karp algorithm.
Specify matching criteria by passing a comma-separated list of column names to ``keys`` and
another comma-separated list of tolerance intervals to ``deltas``. For an exact match,
specify a delta of 0. The order of
``delta`` values must match the order of the ``keys``. Specify the desired matching
ratio via the ``goal`` input; if the matching algorithm cannot achieve the desired
ratio, it will issue a warning and report the achieved ratio.

**Example** Match cases to controls with a 1:3 ratio based on sex (exact match)
and age at diagnosis (match within 1 year)::

		maximizeControls --keys="Sex,AgeAtDx" --deltas="0,1" --goal="3" --input="group.csv"

The default indicator of group membership is the genotype column. However, any
column in the group file may be used provided that it contains only the values [0,1].
To specify a column other than genotype, use the ``condition`` argument.

**Example** Match females (sex=1) to males (sex=0) with a 1:1 ratio based on age at
diagnosis (match within 2 years)::

		maximizeControls --condition="sex" --keys="AgeAtDx" --deltas="2" --goal="1" --input="group.csv"

.. note::
    If there are no suitable matches for some case subjects, **these case subjects may
    be removed**, and will not appear in the output group file. A warning will be issued
    when this occurs with details on how many subjects were lost.

mergeGroups
-----------
Merge 2 or more phenotype/group files.

Optional Arguments [default value]:
 * ``--path``:			        Path to all input files and destination of output files [current directory]
 * ``--phenotypefiles``:		List of phenotype file names, separated by +
 * ``--groupfiles``:			List of group file names, separated by +
 * ``--phenotypeout``:			Name of output file for merged phenotype data (must be specified if phenotypefiles specified)
 * ``--groupout``:				Name of output file for merged group data (must be specified if groupfiles specified)

Output:
 * Group file (``groupout``) containing merged group data
 * Phenotype file (``phenotypeout``) containing merged phenotype data

 The grouping tool allows you to merge two or more phenotype files together, and/or two or
 more group files together. It removes any duplicate records in both file types,
 so that the resulting data files are ready to be run through the pyPheWAS Research Tools.


**Example** Merge 2 ICD9 phenotype files together and 2 group files together::

		generateGroups --phenotypefiles="icd9_one.csv+icd9_two.csv" --groupfiles="group_one.csv+group_two.csv" --phenotypeout="new_icd9.csv" --groupout="new_group.csv"
