from distutils.core import setup

setup(
  name = 'pyPheWAS',
  packages = ['pyPheWAS'], # this must be the same as the name above
  version = '4.0.2',
  description = 'MASI Lab Port of PheWAS into Python',
  author = 'MASI Lab',
  author_email = 'bennett.landman@vanderbilt.edu',
  url = 'https://github.com/MASILab/pyPheWAS', # use the URL to the github repo
  download_url = 'https://github.com/MASILab/pyPheWAS/tarball/0.1', # I'll explain this in a second
  keywords = ['PheWAS', 'ICD-9', 'ICD-10', 'EMR', 'CPT'], # arbitrary keywords
  classifiers = [],
  install_requires=[ 'numpy>=1.16.4',
                     'matplotlib',
                     'scipy>=1.2.1',
                     'pandas>=0.24.2',
                     'statsmodels>=0.10.1',
                     'hopcroftkarp',
                     'tqdm',
                     'pathlib',
                     'biopython',
                     'flask',
					 'flask_cors',
                     ],
  package_data={
    '':['resources/*.csv', 'resources/*.txt','Explorer_GUI/*'],
  },
  scripts=['bin/pyPhewasLookup',
           'bin/pyPhewasModel',
           'bin/pyPhewasPlot',
           'bin/pyPhewasPipeline',
           'bin/pyProwasLookup',
           'bin/pyProwasModel',
           'bin/pyProwasPlot',
           'bin/pyProwasPipeline',
           'bin/censorData',
           'bin/convertEventToAge',
           'bin/createGenotypeFile',
           'bin/createPhenotypeFile',
           'bin/maximizeControls',
           'bin/mergeGroups',
           'bin/PubMedQuery',
           'bin/NoveltyAnalysis',
           'bin/pyPhewasExplorer',
           'bin/pyPhewasExplorerBackEnd'
           ],
)
