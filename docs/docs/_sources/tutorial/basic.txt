Basics
======

Purpose
-------

This tutorial is meant to help get a user started with pyPhewas and troubleshoot problems.


What are Phewas codes?
----------------------

Phewas codes categorize each icd9 into a subsection of a condition. There are over 15,000 ICD-9 codes included in the Phewas categorization. All of these are placed into roughly 2000 Phewas codes. These codes are then generated into Phewas plots and can be analyzed to find associations.


File Format
-----------

Phenotype File
^^^^^^^^^^^^^^

==== ====== ==================
id   icd9   *other covariates*
==== ====== ==================
11   790.29 ...
1    580.8  ...
131  786.59 ...
9999 740.2  ...
==== ====== ==================

This is the file format that is required for the phenotype file. This file is processed by pyPhewas into either

 * For the logarithmic regression, all id-phewas combinations and a 0/1 of whether or not they occurred.
 * For the linear regression, all id-phewas combinations and the count of times that the id corresponded to the given phewas codes

Genotype File
^^^^^^^^^^^^^

===== ==================
id    *other covariates*
===== ==================
1     ...
32    ...
131   ...
200   ...
===== ==================

Depending on what you are using pyPhewas for, the above file format may be all that is necessary. However, more often than not, the file format will look as follows:

===== ======== ==================
id    genotype *other covariates*
===== ======== ==================
1     0        ...
32    0        ...
131   1        ...
200   0        ...
===== ======== ==================

The genotype column and the 0/1 denote the presence or absence of some other condition that the patient may have. The 'genotype' is the default covariate that is passed by the Phewas object into the logarithmic and linear regressions. If you would prefer to use other covariates. They must be specified as outlined in the documentation and below. While using the genotype column is not required, it is highly recommended for the use of Phewas.

.. note:: The order of the columns as shown above is not required, but it does include readability for people opening and reading the files that are input into pyPhewas.