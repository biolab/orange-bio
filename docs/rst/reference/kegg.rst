============================================================
KEGG - Kyoto Encyclopedia of Genes and Genomes (:mod:`kegg`)
============================================================

.. py:currentmodule:: Orange.bio.kegg

.. automodule:: Orange.bio.kegg
   :members:
   :member-order: bysource

DBEntry (:mod:`entry`)
----------------------

The :class:`entry.DBEntry` represents a DBGET databas entry.
The individual KEGG Database interfaces below provide their own
specialization for this base class.
 
.. autoclass:: Orange.bio.kegg.entry.DBEntry
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Databases interface (:mod:`databases`)
-------------------------------------------

.. autoclass:: Orange.bio.kegg.databases.DBDataBase
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.GenomeEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.Genome
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.GeneEntry
   :members:
   :exclude-members:
      alt_names
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.kegg.databases.Genes
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.CompoundEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.kegg.databases.Compound
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.ReactionEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.kegg.databases.Reaction
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.EnzymeEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.kegg.databases.Enzyme
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.PathwayEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.kegg.databases.Pathway
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Pathway (:mod:`pathway`)
-----------------------------

.. autoclass:: Orange.bio.kegg.pathway.Pathway
   :members:
   :exclude-members:
      entrys
   :member-order: bysource
   :show-inheritance:


Utilities
---------

.. autoclass:: Orange.bio.kegg.entry.parser.DBGETEntryParser
   :members:
