#!/usr/bin/env python

import sys
import os

try:
    from setuptools import setup, find_packages
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import setup, find_packages


NAME = 'Orange-Bioinformatics'
DOCUMENTATION_NAME = 'Orange Bioinformatics'

VERSION = '2.6.22'

DESCRIPTION = 'Orange Bioinformatics add-on for Orange data mining software package.'
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'contact@orange.biolab.si'
URL = 'http://orange.biolab.si/download'
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
    'orange3 add-on',
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
)

PACKAGE_DATA = {
}

# Backwards compatibility stub. Should be removed by the 2.7 release.
PY_MODULES = ["_bioinformatics"]

SETUP_REQUIRES = (
    'setuptools',
)

INSTALL_REQUIRES = (
    'Orange3' if sys.version_info > (3, ) else 'Orange',
    'setuptools',
    'numpy',
    'scipy',
    'six',
    'genesis-pyapi>=1.1.2',
    "slumber>=0.7.1",
    "requests>=2.11.1",
    "requests-cache>=0.4.12",
    "serverfiles>=0.2",
    # Dependencies which are problematic to install automatically
    #'openbabel-python', # You get bindings together with the openbabel library and not stand-alone
    #'scipy', # Requires Fortran compiler
    #'matplotlib', # Requires that numpy is installed first
)

if sys.version_info > (3, ):
    INSTALL_REQUIRES = INSTALL_REQUIRES + ("pyqtgraph", "AnyQt")

if sys.version_info < (3, 3):
    INSTALL_REQUIRES = INSTALL_REQUIRES + ("backports.unittest_mock",)

if sys.version_info < (3, 4):
    INSTALL_REQUIRES = INSTALL_REQUIRES + ("singledispatch",)

if sys.version_info < (3, 5):
    INSTALL_REQUIRES = INSTALL_REQUIRES + ("typing",)


EXTRAS_REQUIRE = {
    'MOL_DEPICT': (
        'oasa'
    ),
    'NETWORK': (
        'Orange[NETWORK]'
    ),
}

DEPENDENCY_LINKS = (
#    'http://bkchem.zirael.org/download/bkchem-0.13.0.tar.gz',
#    'http://orange.biolab.si/download/bkchem-0.13.0.tar.gz',
    'http://orange.biolab.si/download/oasa-0.13.1.tar.gz',
    'http://bkchem.zirael.org/download/oasa-0.13.1.tar.gz',
)

ENTRY_POINTS = {
    'orange.addons': (
        'bio = orangecontrib.bio',
    ),
    'orange.widgets': (
        ('Bioinformatics = orangecontrib.bio.widgets',
         'Prototypes = orangecontrib.bio.widgets.prototypes')
        if sys.version_info < (3, ) else
        ('Bioinformatics = orangecontrib.bio.widgets3',)
    ),
    'orange.canvas.help': (
        ('intersphinx = orangecontrib.bio.widgets:intersphinx'
        if sys.version_info < (3,) else
        'html-index = orangecontrib.bio.widgets3:WIDGET_HELP_PATH'),
    )
}

NAMESPACE_PACAKGES = ["orangecontrib", "orangecontrib.bio"]

if __name__ == '__main__':
    setup(
        name = NAME,
        version = VERSION,
        description = DESCRIPTION,
        long_description = LONG_DESCRIPTION,
        author = AUTHOR,
        author_email = AUTHOR_EMAIL,
        url = URL,
        license = LICENSE,
        keywords = KEYWORDS,
        classifiers = CLASSIFIERS,
        packages = PACKAGES,
        package_data = PACKAGE_DATA,
        py_modules = PY_MODULES,
        setup_requires = SETUP_REQUIRES,
        install_requires = INSTALL_REQUIRES,
        extras_require = EXTRAS_REQUIRE,
        dependency_links = DEPENDENCY_LINKS,
        entry_points = ENTRY_POINTS,
        namespace_packages=NAMESPACE_PACAKGES,
        include_package_data = True,
        zip_safe = False,
    )
