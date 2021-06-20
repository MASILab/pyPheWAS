Customizing pyPheWAS
====================

Want to use a home-made PheCode map? Wish that pyPheWAS supported linear regression?
Well, this is the guide for you! pyPheWAS is developed by grad students in
the MASI lab, so it may not contain every feature a user could want. Luckily, where
there's a will (and a modest amount of python knowledge), there's a way.

On this page we cover:

* :ref:`Custom PheCode Maps`
* :ref:`Alternative Regression Methods`
* :ref:`Installing Local Changes`

To get started, open a terminal and clone the pyPheWAS GitHub repo:

.. code-block:: bash

    git clone https://github.com/MASILab/pyPheWAS

After making local changes, make sure you install your modified package by following
the instructions in the :ref:`Installing Local Changes` section.

----------

Custom PheCode Maps
-------------------
**1. Formatting the PheCode Map**

The easiest way to integrate a custom map is to match the custom map's formatting
with pyPheWAS's default map. It is important to keep in mind that separate maps
are required for ICD9 and ICD10 codes. The default ICD9-PheCode map is constructed as
follows (`*`\ required columns).

+----------------------------+--------+------------------------------------------------------+
| Column Name                | Type   | Description                                          |
+============================+========+======================================================+
| ICD9\ :sup:`*`             | str    | ICD version 9 code                                   |
+----------------------------+--------+------------------------------------------------------+
| ICD9 String\ :sup:`*`      | str    | Text description for ICD 9 code                      |
+----------------------------+--------+------------------------------------------------------+
| PheCode\ :sup:`*`          | str    | Phenotype Code                                       |
+----------------------------+--------+------------------------------------------------------+
| Phenotype\ :sup:`*`        | str    | Text description for PheCode                         |
+----------------------------+--------+------------------------------------------------------+
| Excl. Phecodes             | str    | Range of PheCodes to exclude when using this code to |
|                            |        | define case subjects                                 |
+----------------------------+--------+------------------------------------------------------+
| Excl. Phenotypes           | str    | Text description for Excl. Phecodes                  |
+----------------------------+--------+------------------------------------------------------+
| category                   | int    | Numeric category                                     |
+----------------------------+--------+------------------------------------------------------+
| category_string            | str    | Text description for category                        |
+----------------------------+--------+------------------------------------------------------+

Notes
	* Only four columns are required; these are sufficient for running the full
	  pyPheWAS pipeline (i.e. Lookup, Model, Plot).
	* Leaving out the optional columns *Excl. Phecodes* and *Excl. Phenotypes* will
	  disable the ``--exclude_phecode=map`` option in the :ref:`createPhenotypeFile` tool.
	* Leaving out the optional columns *category* and *category_string* will disable
	  category coloring in the Manhattan and Log Odds plots produced by :ref:`pyPhewasPlot`.
	* The ICD10-PheCode map is identical to above, with the exception
	  of columns *ICD9* and *ICD9 String* being replaced by *ICD10* and *ICD10 String*,
	  respectively.


**2. Incorporating the map into pyPheWAS**

The only necessary modification in this case will be telling pyPheWAS where to find the map file.
First, **save your custom PheCode map in the resources folder** of your local clone of the
pyPheWAS repository: ``pyPheWAS/pyPheWAS/resources/``. This is also where the default maps
are stored.

Next, tell pyPheWAS to **load your custom map as either the ICD9, ICD10, or CPT phecode map**.
This is done at the end of the `pyPhewasCorev2.py:L1154`_ script, as shown below.


.. code-block:: python

    #-----------------------------------------------------------------------------
    # load ICD maps (pyPheWAS)
    # Change the filename argument to match the name of your custom PheCode map(s)
    icd9_codes = get_codes('phecode_map_v1_2_icd9.csv')
    icd10_codes = get_codes('phecode_map_v1_2_icd10_beta.csv')
    # load CPT maps (pyProWAS)
    # Change the filename argument to match the name of your custom ProCode map
    cpt_codes = get_codes('prowas_codes.csv')
    #-----------------------------------------------------------------------------

Notes
    * The ICD9 and ICD10 PheCode maps are merged to obtain one master PheCode list.
    * To install your changes and continue using the command line pyPheWAS
      interface, make sure you follow the instructions in the
      :ref:`Installing Local Changes` section.


Alternative Regression Methods
------------------------------
Changes to the regression method may be implemented by modifying the function
:py:func:`pyPheWAS.pyPhewasCorev2.fit_pheno_model`. This function accepts feature
arrays and other regression settings and returns a statistics vector for the
fitted regression. Though technically any statistics package may be used to implement
an alternative regression method, we recommend using the
`statsmodels <https://www.statsmodels.org/stable/index.html>`_
package.

Once you have picked the alternative model from statsmodels, you will need to make
several modifications to ``fit_pheno_model``.

.. note::

    The ``fit_pheno_model`` function contains two distinct regression fits:
    one with regularization (``lr=1``) and one without regularization (``lr=0``).
    The implementation details of these methods vary slightly. However,
    the fit without regularization is only used when the ``--legacy`` flag is
    activated in :ref:`pyPhewasModel`, so *only the regularized fit implementation
    will be described here.*

**1. Setting up the regression function**

The Logit regression function object is declared by passing the constructor
an array of response variable values and a matrix of predictor values (as
shown in the snippet from `pyPhewasCorev2.py:L412`_ below). Modify this line to
match the declaration structure of your alternate regression.

.. code-block:: python

    # column 'y' is the PheCode vector
    predictors = covariates.replace(" ", "").split('+')
    predictors[0] = 'y'
    f = [response.strip(), predictors]
    logit = sm.Logit(data[f[0]], data[f[1]])


**2. Fitting the regression function**

The next line fits the regression function with regularization (`pyPhewasCorev2.py:L413`_).
Again, modify this as needed to match your alternate regression method.

.. code-block:: python

    model = logit.fit_regularized(method='l1', alpha=0.1, disp=0, trim_mode='size', qc_verbose=0)


**3. Formatting the regression function stats**

Finally, you need to pull the stats from your fitted regression. ``fit_pheno_model``
should return the following values (in order):
-log\ :sub:`10`\ (p-value),
p-value, beta, beta’s confidence interval, and beta’s standard error.
These are pulled from the fitted model as shown below (`pyPhewasCorev2.py:L417`_).
Check the API of your alternative regression model to ensure that these
values are the same.

.. code-block:: python

    # get results for y (the PheCode vector)
    p = model.pvalues.y
    beta = model.params.y
    conf = model.conf_int()
    conf_int = '[%s,%s]' % (conf[0]['y'], conf[1]['y'])
    stderr = model.bse.y
    reg_result = [-math.log10(p), p, beta, conf_int, stderr]  # collect results


Installing Local Changes
------------------------
After making changes to the pyPheWAS repository, you may install your local version
of the package by running the following from a terminal:


.. code-block:: bash

    cd pyPheWAS # change to the top-level pyPheWAS repository
    python setup.py sdist # build the local package
    pip install . --upgrade # install the package



.. _pyPhewasCorev2.py:L1154: https://github.com/MASILab/pyPheWAS/blob/f7cae0756dd2792cdf9bf166446c8d75ed33d972/pyPheWAS/pyPhewasCorev2.py#L1154
.. _pyPhewasCorev2.py:L412: https://github.com/MASILab/pyPheWAS/blob/f7cae0756dd2792cdf9bf166446c8d75ed33d972/pyPheWAS/pyPhewasCorev2.py#L412
.. _pyPhewasCorev2.py:L413: https://github.com/MASILab/pyPheWAS/blob/f7cae0756dd2792cdf9bf166446c8d75ed33d972/pyPheWAS/pyPhewasCorev2.py#L413
.. _pyPhewasCorev2.py:L417: https://github.com/MASILab/pyPheWAS/blob/f7cae0756dd2792cdf9bf166446c8d75ed33d972/pyPheWAS/pyPhewasCorev2.py#L417
