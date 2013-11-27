.. py:currentmodule:: Orange.bio.obiDictyMutants

.. index:: mutant phenotypes
.. index:: phenotypes
.. index::
   single: D. dictyostelium; mutants

********************************************************
D. discoideum Mutant Phenotypes (:mod:`obiDictyMutants`)
********************************************************

Module :mod:`obiDictyMutants` provides interface to `Dictyostelium mutant phenotypes
<http://dictybase.org/Downloads/>`_ data from the dictyBase_.
The mutants are presented as `DictyMutant` objects with their respective name,
strain descriptor, associated genes and associated phenotypes.

Classes and Functions
=====================

.. automodule:: Orange.bio.obiDictyMutants
   :members:
   :member-order: bysource

Examples
========

   >>> from Orange.bio.obiDictyMutants import *
   >>> # Create a set of all mutant objects
   >>> dicty_mutants = mutants()
   >>> # List a set of all genes referenced by a single mutant
   >>> print mutant_genes(dicty_mutants[0])
   ['cbfA']
   >>> # List a set of all phenotypes referenced by a single mutant
   >>> print mutant_phenotypes(dicty_mutants[0])
   ['aberrant protein localization']
   >>> # List all genes or all phenotypes referenced on Dictybase
   >>> print genes()
   >>> print phenotypes()
   >>> # Display a dictionary {phenotypes: set(mutant_objects)}
   >>> print phenotype_mutants()
   >>> # Display a dictionary {genes: set(mutant_objects)}
   >>> print gene_mutants()

