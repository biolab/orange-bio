#!/usr/bin/env python
"""\
Orange Bioinformatics
=====================

Orange Bioinformatics is an add-on for Orange data mining 
software package. It extends Orange by providing common functionality
for basic tasks in bioinformatics.
"""

DOCLINES = __doc__.split("\n")

try:
    from setuptools import setup
    have_setuptools = True
except ImportError:
    from distutils.core import setup
    have_setuptools = False

import os, glob

CLASSIFIERS = """\
Development Status :: 4 - Beta
Programming Language :: Python
License :: OSI Approved :: GNU General Public License (GPL)
Operating System :: POSIX
Operating System :: Microsoft :: Windows
Topic :: Scientific/Engineering :: Visualization
Topic :: Scientific/Engineering :: Bio-Informatics
Intended Audience :: Education
Intended Audience :: Science/Research
"""

KEYWORDS = """\
data mining 
machine learning,
artificial intelligence
bioinformatics,
gene ontology
KEGG
expression profiles
"""                      

NAME                = "Orange-Bioinformatics"
DESCRIPTION         = DOCLINES[0]
LONG_DESCRIPTION    = "\n".join(DOCLINES[3:])
URL                 = "http://www.biolab.si/obi/"
DOWNLOAD_URL        = "https://bitbucket.org/biolab/orange-addon-bioinformatics/downloads"
LICENSE             = "GNU General Public License (GPL)"
CLASSIFIERS         = filter(None, CLASSIFIERS.split("\n"))
AUTHOR              = "Bioinformatics Laboratory, FRI UL"
AUTHOR_EMAIL        = "orange@fri.uni-lj.si"
KEYWORDS            = filter(None, KEYWORDS.split('\n'))

MAYOR = 1
MINOR = 1
MICRO = 0
ISRELEASED = False

VERSION = "{0}.{1}a.{2}".format(MAYOR,MINOR, MICRO)


# list all documentation files that need to be included
docFiles = []
for (dirp, dirns, n) in os.walk('doc'):
    nr = [n1.replace('\\', '/') for n1 in n]
    dirn = dirp.replace('\\', '/')[4:]
    if len(dirn):
        dirn = dirn + '/'
    docFiles.extend( [dirn + n1r for n1r in nr if '.svn' not in dirp + '/' + n1r] )

DEST_DIR="Orange/add-ons/Bioinformatics"


if os.path.exists("VERSION.txt"):
    VERSION = open("VERSION.txt", "rb").read()

if have_setuptools:
    setuptool_args = {"install_requires": ["Orange", "suds"],
                      "zip_safe": False,
                     }
else:
    setuptool_args = {}
    
PACKAGES = [ 'widgets', 'widgets.prototypes', 'doc', '',
             'obiKEGG2', 'obiKEGG2.entry' ]
    
PACKAGE_DATA = {'widgets': ['icons/*.png'],
                'doc': docFiles,
                '':["addon.xml"]
                }
                          
if __name__ == "__main__":
    setup(name = NAME,
          version = VERSION,
          description = DESCRIPTION,
          long_description = LONG_DESCRIPTION,
          author = AUTHOR,
          author_email = AUTHOR_EMAIL,
          url = URL,
          download_url = DOWNLOAD_URL,
          license = LICENSE,
          keywords = KEYWORDS,
          classifiers = CLASSIFIERS,
          package_dir = {"": "."},
          packages = PACKAGES,
          package_data = PACKAGE_DATA,
          extra_path=("orange-bioinformatics", DEST_DIR),
          
          **setuptool_args)
