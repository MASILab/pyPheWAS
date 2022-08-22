## PyPheWAS

Repository for the pyPheWAS project.
Full documentation at https://pyphewas.readthedocs.io/en/latest/

### Developers
Cailey Kerley, PhD Candidate

Shikha Chaganti, PhD

Bennett Landman, PhD

## Cite pyPheWAS
Kerley, C.I., Chaganti, S., Nguyen, T.Q. et al. pyPheWAS: A Phenome-Disease Association Tool for Electronic Medical Record Analysis. *Neuroinform* (2022). https://doi.org/10.1007/s12021-021-09553-4


## Latest Release: pyPheWAS 4.1

#### 4.1.0
- pyPheWAS Explorer updates
- New demographic variables added to synthetic dataset

#### 4.0.5
- convertEventToAge includes new warning for calculated ages are negative
- small bugs fixed in maximizeControls, NoveltyAnalysis, and PubMedQuery tools


#### 4.0.4
- createPhenotypeFile now supports more options for controlling case/control group curation
- Documentation updates

#### 4.0.3
- **Novelty Analysis** tools: examine the relative literary novelty of disease-phecode pairings
- **pyPheWAS Explorer**: an interactive visualization of PheDAS experiments
- createGenotypeFile updated - now called createPhenotypeFile
- Minor bug fixes


### Older Releases

#### pyPheWAS 3.2.0
- Configurable threshold for number of subjects required to run the regression on an individual PheCode
- All regressions are now fit with regularization (old scheme available with 'legacy' option)
- Minor changes to Manhattan plot
- PheCode/ProCode categories added to regression file
- Minor bug fixes

#### pyPheWAS 3.1.1
- New Analysis Type: ProWAS Tools
- New Plot Type: Volcano Plot (see pyPhewasPlot)
- maximizeControls now saves explicit Case/Control matches
- New PheCode category colors in plots are more distinguishable
- Improved command line tool argument handling
- Improved error handling
- Documentation overhaul
- API update
- Minor bug fixes

#### pyPheWAS 3.0.1
- Bug fixes including __FDR & Bonferroni threshold calculations__
- Header saved in feature matrices
- More file formats available for saving plots

#### pyPheWAS 3.0.0
- Support for both ICD 9 and ICD 10
- All 3 regression types (binary, count, & duration) optimized for big data
- pyPhewasPipeline: a streamlined combination of pyPhewasLookup, pyPhewasModel, and pyPhewasPlot
- Compatibility with Python 3
- Age matching now saves the explicit mapping between controls/cases in addition to the resulting group file
- Operation of the ICD censoring function matches the description in the documentation
- Minor bug fixes
