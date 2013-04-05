============================================================
KEGG - Kyoto Encyclopedia of Genes and Genomes (:mod:`kegg`)
============================================================


.. automodule:: Orange.bio.obiKEGG
   :members:
   :member-order: bysource

DBEntry (:mod:`entry`)
----------------------

The :class:`~.entry.DBEntry` represents a DBGET databas entry.
The individual KEGG Database interfaces below provide their own
specialization for this base class.
 
.. autoclass:: Orange.bio.obiKEGG.entry.DBEntry
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Databases interface (:mod:`databases`)
-------------------------------------------

.. autoclass:: Orange.bio.obiKEGG.databases.DBDataBase
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.GenomeEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.Genome
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.GeneEntry
   :members:
   :exclude-members:
      alt_names
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.obiKEGG.databases.Genes
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.CompoundEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.obiKEGG.databases.Compound
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.ReactionEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.obiKEGG.databases.Reaction
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.EnzymeEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: Orange.bio.obiKEGG.databases.Enzyme
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.PathwayEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: Orange.bio.obiKEGG.databases.Pathway
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Pathway (:mod:`pathway`)
-----------------------------

.. autoclass:: Orange.bio.obiKEGG.pathway.Pathway
   :members:
   :exclude-members:
      entrys
   :member-order: bysource
   :show-inheritance:


Utilities
---------

.. autoclass:: Orange.bio.obiKEGG.entry.parser.DBGETEntryParser
   :members:
