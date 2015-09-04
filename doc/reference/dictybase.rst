.. py:currentmodule:: Orange.bio.dicty.phenotypes

.. index:: mutant phenotypes
.. index:: phenotypes
.. index::
   single: D. dictyostelium; mutants

*********************************************************
D. discoideum Mutant Phenotypes (:mod:`dicty.phenotypes`)
*********************************************************

This modules provides an interface to `Dictyostelium mutant
phenotypes <http://dictybase.org/Downloads/>`_ data from the
`dictyBase <http://dictybase.org/>`_.  The mutants are presented as
:obj:`DictyMutant` objects with their respective name, strain descriptor,
associated genes and associated phenotypes.

>>> from Orange.bio.dicty.phenotypes import *
>>> # Create a set of all mutant objects
>>> dicty_mutants = mutants()
>>> # List a set of all genes referenced by a single mutant
>>> print mutant_genes(dicty_mutants[0])
['cbfA']
>>> # List a set of all phenotypes referenced by a single mutant
>>> print mutant_phenotypes(dicty_mutants[0])
['aberrant protein localization']

Classes and Functions
=====================

.. automodule:: Orange.bio.dicty.phenotypes
   :members:
   :member-order: bysource


