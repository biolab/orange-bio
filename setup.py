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

icons = glob.glob(os.path.join('widgets', 'icons', '*.png'))

setup(name = "Genomics",
      version = "VERSION",
      description = "Genomics extensions for Orange",
      author="University of Ljubljana, AI lab",
      author_email="tomaz.curk@fri.uni-lj.si",
      packages = [ 'widgets' ],
      data_files = [('doc', docFiles),
                    (os.path.join('widgets', 'icons'), icons)],
      extra_path="Genomics",
      py_modules = [ 'orngKEGG', 'orngGsea', 'orngGeneMatcher' ],
      scripts=["registerWidgets.py"]
      )
