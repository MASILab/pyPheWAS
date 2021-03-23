pyPheWAS Explorer Setup
=======================

Exactly like the :ref:`PheWAS Tools`, the Explorer
requires a **group file** and a **phenotype file**, the formats of which
may be found in the :ref:`File Formats` section. These files must be called
``group.csv`` and ``icds.csv``, respectively, and be saved within the same directory
on your filesystem.

.. note:: At this time, the Explorer can only parse numeric group variables. Categorical
   variables (e.g. sex, race) must be converted to numeric categories before launching
   pyPheWAS Explorer.

Launching the Explorer
----------------------
The script ``pyPhewasExplorer`` is installed with the pyPheWAS package to launch
the Explorer. This script accepts only two input arguments:

* ``--indir``:		Path to input directory [default: current directory]
* ``--response``:	Column name from group.csv to use as the dependent variable [default: target]

**Example**::

		pyPhewasExplorer --indir /Users/me/Documents/EMRdata/ --response ADHD

To launch the Explorer, run the launch script as shown above and then open
``http://localhost:8000/`` in a web browser (preferably Google Chrome). 
