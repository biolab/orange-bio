from distutils.core import setup
import os, glob

# list all documentation files that need to be included
docFiles = []
for (dirp, dirns, n) in os.walk('doc'):
    nr = [n1.replace('\\', '/') for n1 in n]
    dirn = dirp.replace('\\', '/')[4:]
    if len(dirn):
        dirn = dirn + '/'
    docFiles.extend( [dirn + n1r for n1r in nr if '.svn' not in dirp + '/' + n1r] )

destDir="Orange/add-ons/Bioinformatics"

if __name__ == "__main__":
    setup(name = "Orange-Bioinformatics",
          version = "1.0.0b",
          description = "Bioinformatics Add-On for Orange",
          author = "Bioinformatics Laboratory, FRI UL",
          author_email = "orange@fri.uni-lj.si",
          url = "http://www.biolab.si/obi/",
          download_url = "https://bitbucket.org/biolab/orange-addon-bioinformatics",
          packages = [ 'widgets', 'widgets.prototypes', 'doc', '.' ],
          package_data = {'widgets': ['icons/*.png'], 'doc': docFiles, '.':["addon.xml"] },
          extra_path=("orange-bioinformatics", destDir),
          license = "GNU General Public License (GPL)",
          keywords = ["data mining", "machine learning", "artificial intelligence",
                      "bioinformatics", "gene ontology", "KEGG", "expression profiles"],
          classifiers = ["Development Status :: 4 - Beta",
                     "Programming Language :: Python",
                     "License :: OSI Approved :: GNU General Public License (GPL)",
                     "Operating System :: POSIX",
                     "Operating System :: Microsoft :: Windows",
                     "Topic :: Scientific/Engineering :: Visualization",
                     "Topic :: Scientific/Engineering :: Bio-Informatics",
                     "Intended Audience :: Education",
                     "Intended Audience :: Science/Research"
                     ],
          long_description="""\
Orange Bioinformatics
=====================

Orange Bioinformatics is an add-on for Orange data mining 
software package. It extends Orange by providing common functionality
for basic tasks in bioinformatics.

""")
