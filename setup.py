#!/usr/bin/env python

try:
    import distribute_setup
    distribute_setup.use_setuptools()
except ImportError:
    # For documentation we load setup.py to get version
    # so it does not matter if importing fails
    pass

import os

from setuptools import setup, find_packages

NAME = 'Orange-Bioinformatics'
DOCUMENTATION_NAME = 'Orange Bioinformatics'

VERSION = '2.5a4'

DESCRIPTION = 'Orange Bioinformatics add-on for Orange data mining software package.'
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'contact@orange.biolab.si'
URL = 'http://orange.biolab.si/addons/'
DOWNLOAD_URL = 'https://bitbucket.org/biolab/orange-bioinformatics/downloads'
DOCS_URL = 'http://orange-bioinformatics.readthedocs.org/'
LICENSE = 'GPLv3'

KEYWORDS = (
    'data mining',
    'machine learning',
    'artificial intelligence',
    'bioinformatics',
    'gene ontology',
    'KEGG',
    'expression profiles',
    'microarray',
    'genomics',
    'orange',
    'orange add-on',
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

SETUP_REQUIRES = (
    'setuptools',
)

INSTALL_REQUIRES = (
    'Orange',
    'setuptools',
    'numpy',
    # Dependencies which are problematic to install automatically
    #'openbabel-python', # You get bindings together with the openbabel library and not stand-alone
    #'scipy', # Requires Fortran compiler
    #'matplotlib', # Requires that numpy is installed first
),

EXTRAS_REQUIRE = {
    'GUI': (
        # Dependencies which are problematic to install automatically
        #'PyQt', # No setup.py
    ),
    'MOL_DEPICT': (
        'oasa'
    ),
    'NETWORK': (
        'Orange[NETWORK]'
    ),
    'KEGG': (
        'suds'
    )

}

DEPENDENCY_LINKS = (
#    'http://bkchem.zirael.org/download/bkchem-0.13.0.tar.gz',
#    'http://orange.biolab.si/download/bkchem-0.13.0.tar.gz',
    'http://bkchem.zirael.org/download/oasa-0.13.1.tar.gz',
    'http://orange.biolab.si/download/oasa-0.13.1.tar.gz',
)

ENTRY_POINTS = {
    'orange.addons': (
        'bioinformatics = _bioinformatics',
    ),
    'orange.widgets': (
        'Bioinformatics = _bioinformatics.widgets',
        # This should be unneeded, because module given should load (register)
        # all wanted widgets and prototypes should just have a flag, but for now ...
        'Prototypes = _bioinformatics.widgets.prototypes',
    ),
}

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
        docs_url = DOCS_URL,
        license = LICENSE,
        keywords = KEYWORDS,
        classifiers = CLASSIFIERS,
        packages = PACKAGES,
        package_data = PACKAGE_DATA,
        setup_requires = SETUP_REQUIRES,
        install_requires = INSTALL_REQUIRES,
        extras_require = EXTRAS_REQUIRE,
        dependency_links = DEPENDENCY_LINKS,
        entry_points = ENTRY_POINTS,
        include_package_data = True,
        zip_safe = False,
    )
