from distutils.core import setup
setup(
  name = 'pyPheWAS',
  packages = ['pyPheWAS'], # this must be the same as the name above
  version = '0.1.1',
  description = 'MASI Lab Port of PheWAS into Python',
  author = 'MASI Lab',
  author_email = 'bennett.landman@vanderbilt.edu',
  url = 'https://github.com/bennettlandman/pyPheWAS', # use the URL to the github repo
  download_url = 'https://github.com/bennettlandman/pyPheWAS/tarball/0.1', # I'll explain this in a second
  keywords = ['PheWAS', 'ICD-9', 'EMR'], # arbitrary keywords
  classifiers = [],
  install_requires=['numpy',
	'matplotlib',
	'scipy',
	'pandas',
	'statsmodels',
	],
)
