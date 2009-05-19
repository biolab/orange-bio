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

destDir="orange/add-ons/Bioinformatics"

if __name__ == "__main__":
    setup(name = "Orange Bioinformatics",
          version = "1.0",
          description = "Bioinformatics Add-On for Orange",
          author="University of Ljubljana, AI lab",
          maintainer_email="tomaz.curk@fri.uni-lj.si",
          packages = [ 'widgets', 'widgets.prototypes', 'doc' ],
          package_data = {'widgets': ['icons/*.png'], 'doc': docFiles},
          extra_path=("orange-bioinformatics", destDir),
          py_modules = ['obiKEGG', 'obiGsea', 'obiData', 'obiGenomicsUpdate', 'stats', 'pstat', 'obiExpression', 'obiGO', 'obiProb', 'obiAssess', 'obiGeneSets', 'obiMeSH', 'obiDicty', 'obiTaxonomy', 'obiChem', 'obiGene', 'obiGEO' ],
          scripts=["post_install_script.py"]
          )
