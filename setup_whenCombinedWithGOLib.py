from distutils.core import setup, Extension
import os, glob

# delete files on exclude list
for fn in open("exclude.lst", "rt").readlines():
	try:
		os.remove(fn.rstrip('\n\r'))
	except:
		pass

# list all documentation files that need to be included
docFiles = []
for (dirp, dirns, n) in os.walk('doc'):
	nr = [n1.replace('\\', '/') for n1 in n]
	dirn = dirp.replace('\\', '/')[4:]
	if len(dirn):
		dirn = dirn + '/'
	docFiles.extend( [dirn + n1r for n1r in nr if '.svn' not in dirp + '/' + n1r] )

module=Extension("_GOLib", sources=["go.c"])

setup(name = "Genomics",
      version = "0.9.67",
      description = "Genomics extensions for Orange",
      author="University of Ljubljana, AI lab",
      author_email="tomaz.curk@fri.uni-lj.si",
      ext_modules=[module],
      packages = [ 'widgets', 'doc', 'data' ],
      package_data = {'widgets': ['icons/*.png'], 'doc': docFiles},
      extra_path="Genomics",
      py_modules = [ 'go', 'obiKEGG', 'obiGsea', 'obiGeneMatch', 'obiData', 'obiGenomicsUpdate', 'obiExpression', 'stats', 'pstat' ],
      scripts=["registerWidgets.py", "post_install_script.py"]
      )
