Basics
======

Purpose
-------

This tutorial is meant to help get a user started with pyPhewas and troubleshoot problems.

Installation
------------

pyPheWAS is compatible with Python 2.7+ and Python 3. Install the pyPheWAS package by running::

		pip install pyPheWAS

If this command fails, make sure that you are running Python 2.7+ and that you have pip set up on your machine.

As long as the install is successful, the pyPheWAS package can now be run from any directory.

.. note:: If installing on a computing cluster (or other environment in which you do not have admin privileges) it may necessary to install pyPheWAS locally using pip's *--user* flag.


What are Phewas codes?
----------------------

**PheWAS**:  Phenome Wide Association Studies

PheWAS codes categorize ICD codes into groups of related codes. There are over 15,000 ICD-9 and ICD-10 codes
included in the Phewas categorization. All of these are mapped to 1865 PheWAS codes. The mappings used for this package
and more information about PheWAS codes can be found `here <https://phewascatalog.org/>`_. The ICD-9 map currently used
used is `version 1.2 <https://phewascatalog.org/phecodes>`_ and the ICD-10 map is
`version 1.2b1 <https://phewascatalog.org/phecodes_icd10>`_.


File Formats
------------

Phenotype File
^^^^^^^^^^^^^^

==== ======== ======== ========
id   ICD_CODE ICD_TYPE AgeAtICD
==== ======== ======== ========
11   790.29   9        10.4
11   580.8    9        11.5
131  A03.2    10       60.0
9999 740.2    9        0.2
==== ======== ======== ========

This is the file format that is required for the phenotype file. This file is processed by pyPhewas into either

 * For the logarithmic regression, all id-phewas combinations and a 0/1 of whether or not they occurred.
 * For the linear regression, all id-phewas combinations and the count of times that the id corresponded to the given phewas codes

Genotype File
^^^^^^^^^^^^^

===== ======== ==================
id    genotype *other covariates*
===== ======== ==================
1     0        ...
32    0        ...
131   1        ...
200   0        ...
===== ======== ==================

The genotype column and the 0/1 denote the presence or absence of some other condition that the patient may have. The
'genotype' is the default prediction variable that is passed by the Phewas object into the logarithmic and linear
regressions. If you would prefer to use other covariates, they must be specified as outlined in the documentation. While
using the genotype column is not required, it is highly recommended for the use of Phewas.

.. note:: The order of the columns as shown above is not required, but it does encourage readability when opening and reading the files that are input into pyPhewas.
