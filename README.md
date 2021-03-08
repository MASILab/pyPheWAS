## PyPheWAS

Repository for the pyPheWAS project.
Full documentation at https://pyphewas.readthedocs.io/en/latest/

### Developers
Cailey Kerley, B.S.

Shikha Chaganti, PhD

Bennett Landman, PhD

## Latest Release: pyPheWAS 4.0.1
This release includes:
- **Novelty Analysis** tools: examine the relative literary novelty of disease-phecode pairings
- **pyPheWAS Explorer**: an interactive visualization of PheDAS experiments
- createGenotypeFile updated - now called createPhenotypeFile
- Minor bug fixes


### Older Releases

### pyPheWAS 3.2.0
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
