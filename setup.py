#!/usr/bin/env python

import os

from setuptools import setup, find_packages

NAME = 'Orange-Bioinformatics'

VERSION = '1.1a'

DESCRIPTION = 'Orange Bioinformatics is an add-on for Orange data mining software package.'
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'contact@orange.biolab.si'
URL = 'http://orange.biolab.si/addons/'
DOWNLOAD_URL = 'https://bitbucket.org/biolab/orange-bioinformatics/downloads'
LICENSE = 'GPLv3'

KEYWORDS = (
    'data mining',
    'machine learning',
    'artificial intelligence',
    'bioinformatics',
    'gene ontology',
    'KEGG',
    'expression profiles',
)

CLASSIFIERS = (
    'Development Status :: 4 - Beta',
    'Environment :: X11 Applications :: Qt',
    'Environment :: Console',
    'Environment :: Plugins',
    'Programming Language :: Python',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'Topic :: Scientific/Engineering :: Visualization',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
)

PACKAGES = find_packages(
    exclude = ('*.tests', '*.tests.*', 'tests.*', 'tests'),
)
 
PACKAGE_DATA = {
}

INSTALL_REQUIRES = (
    'Orange',
    'suds',
    'numpy',
    'requests',
    'scipy',
    'oasa',
    'bkchem',
    'matplotlib',
    'PIL',
    'sqlite3',
    'networkx',
    'pygraphviz',
    'PyQt4',
    'openbabel',
),

DEPENDENCY_LINKS = (
    'http://bkchem.zirael.org/download/bkchem-0.13.0.tar.gz',
    'http://bkchem.zirael.org/download/oasa-0.13.1.tar.gz',
)

ENTRY_POINTS = {
    'orange.widgets': (
        'bioinformatics = Orange.bioinformatics.widgets',
    ),
}

NAMESPACE_PACKAGES = (
    'Orange',
)

if __name__ == '__main__':
    setup(
        name = NAME,
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
        packages = PACKAGES,
        package_data = PACKAGE_DATA,
        install_requires = INSTALL_REQUIRES,
        dependency_links = DEPENDENCY_LINKS,
        entry_points = ENTRY_POINTS,
        namespace_packages = NAMESPACE_PACKAGES,
        include_package_data = True,
        zip_safe = False,
    )
