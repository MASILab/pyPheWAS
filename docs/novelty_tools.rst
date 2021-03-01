Novelty Analysis Tools
======================
This page describes the command line tools available for running a Novelty
Analysis.

These tools require a **regression** file, which is described in the
:ref:`pyPhewasModel` section.

Overview of Novelty Analysis Tools
-----------------------------------
pyPheWAS's Novelty Analysis tools are a post-processing step for PheDAS analyses
that aim to estimate the relative "novelty" of a disease-phecode association. In brief,
this involves calculating a **Novelty Finding Index (NFI)** for each disease-phenotype association
that measures the *degree to which it is already known* based on data mined from PubMed abstracts.
If a disease-phenotype pairing is present in a large number of PubMed abstracts, the
association is assigned a low NFI and considered well known.
In contrast, if a disease-phenotype pairing is present in only a few PubMed abstracts,
the association is assigned a high NFI and considered relatively unknown.
For more information on the NFI, please see our publication [Chaganti2019b].

Novelty Analysis Functions:

* :ref:`PubMedQuery`: run a PubMed search for A) all PheCodes or B) a set
  of custom terms (i.e. a disease)
* :ref:`NoveltyAnalysis`: calculate and visualize the Novelty Finding Index for
  the results from a PheDAS


PubMedQuery
-----------
Run a PubMed Search.

Conducts a search of PubMed Titles, Abstracts, and Keywords for A) all 1866 PheCodes
in the pyPheWAS map, or B) a list of custom search terms. Search results are saved to a
CSV file, where the 'IdsList' column includes a list of unique PubMed article identifiers
found for each search. Both search types must be performed in order to proceed to the
:ref:`NoveltyAnalysis` step.

*Option A: Mass PheCode Search*

This option iterates over all 1,866 PheCodes in the PheWAS mapping and searched PubMed
for related articles. (It takes approximately 24 hours to conduct searches for all 1,866 PheCodes.)
The mass search is done by mapping a PheCode's corresponing ICD-9
and ICD-10 codes to CUIs (Concept Unique Identifiers) in the
`UMLS Metathesaurus <https://www.nlm.nih.gov/research/umls/knowledge_sources/metathesaurus/index.html>`.
Then, all strings that correspond to those CUIs are used to query PubMed.
**Warning**: due to the size of the UMLS Metathesaurus, this option requires a
machine with at least 16 GB of RAM.

.. note::
   A free license is required to download the UMLS Metathesaurus. For more information, please
   see the `UMLS website <https://www.nlm.nih.gov/research/umls/index.html>`.

*Option B: Custom Search*

This option searches the PubMed database for a set of custom search strings, specified
via a text file with one search string per line. (See the examples below for an example
file.)

Required Arguments:
 * ``--outdir``: 	Path to output directory

Optional Arguments (at least one of these must be specified):
 * ``--umls``:		     Path to the UMLS Metathesaurus file (file should be called 'MRCONSO.RRF')
 * ``--custom_terms``: Path to a file containing custom search terms (i.e for the target disease)


Output:
 CSV file(s) that contain a list of PubMed article IDs for each query (either a PheCode
 or a custom search).

**Example** Run a mass PheCode PubMed search::

		PubMedQuery --outdir="/Users/me/Documents/PubMed_data/" --umls="/Users/me/Documents/MRCONSO.RRF"

**Example** Run a custom PubMed search on terms in **ADHD_terms.txt**.

**ADHD_terms.txt**::

  | ADHD |
  | ADDH |
  | attention deficit hyperactivity disorder |
  | attention deficit disorder with hyperactivity |
  | attention deficit |

**Command**::

		PubMedQuery --outdir="/Users/me/Documents/EMRdata" --custom_terms="/Users/me/Documents/ADHD_terms.txt"


NoveltyAnalysis
---------------
Calculate the Novelty Finding Index for a set of PheDAS results.

Calculates the NFI for all PheCodes in a regression file output from :reg:`pyPhewasModel`.
All calculated values are saved to the regression, and novelty plots are created
for all PheCodes with a second generation p-value of 0. For more information on how
the NFI is calculated, please see our publication [Chaganti2019b].

Required Arguments:
 * ``--pm_dir``:     Path to directory where mass PheCode PubMed search results are stored
 * ``--statfile``:   Name of the regression file output from :reg:`pyPhewasModel`
 * ``--dx_pm``:		   Name of the disease's PubMed search results file (obtained via :ref:`PubMedQuery` custom search)
 * ``--null_int``:   Null interval to use in calculating the NFI

Optional Arguments [default value]:
 * ``--path``:			Path to all input files and destination of output files [current directory]
 * ``--postfix``:	  Descriptive postfix for output files (e.g. poster or ages50-60)


Output:
 NFI calculations saved with the regression file and novelty plots for significant
 (2nd generation pvalue = 0) PheCodes.

 Additional regression file columns include:
  * sgpv: second generation p-value
  * ppv: positive predictive value
  * ecdf: empirical cumulative distribution function estimated from the PubMed Proportions
  * DX_PM_count: number of PubMed results found for the target disease
  * phe_PM_count: number of PubMed results found for each PheCode
  * joint_PM_count: number of PubMed results that mention both the target disease and a PheCode
  * P_PM_phe_given_DX: PubMed Proportion
  * Novelty_Finding_Index: the NFI for each PheCode

**Example** Calculate the NFI for a PheDAS regression of ADHD::

		NoveltyAnalysis --null_int="[0.3,1.1]" --pm_dir="/Users/me/Documents/PubMed_data/" --dx_pm="ADHD_pubmed_search.csv" --statfile="regressions.csv" --path="/Users/me/Documents/EMRdata/"


.. note::
   The null interval (`null_int`) is specified in terms of the odds ratio, but
   results are plotted using the log odds ratio.
