from distutils.core import setup, Extension
import os, glob

# list all documentation files that need to be included
docFiles = []
for (dirp, dirns, n) in os.walk('doc'):
	nr = [n1.replace('\\', '/') for n1 in n]
	dirn = dirp.replace('\\', '/')[4:]
	if len(dirn):
		dirn = dirn + '/'
	docFiles.extend( [dirn + n1r for n1r in nr if '.svn' not in dirp + '/' + n1r] )

setup(name = "Genomics",
      version = "1.0",
      description = "Genomics extensions for Orange",
      author="University of Ljubljana, AI lab",
      author_email="tomaz.curk@fri.uni-lj.si",
      packages = [ 'widgets', 'doc' ],
      package_data = {'widgets': ['icons/*.png'], 'doc': docFiles},
      extra_path="Genomics",
      py_modules = [ 'obiKEGG', 'obiGsea', 'obiGeneMatch', 'obiData', 'obiGenomicsUpdate', 'obiExpression', 'stats', 'pstat' ],
      scripts=["registerWidgets.py"]
      )
