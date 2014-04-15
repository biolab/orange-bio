============================================================
Protein-protein interactions (:mod:`ppi`)
============================================================

.. py:currentmodule:: Orange.bio.ppi

:class:`PPIDatabase` is an abstract class defining a common interface
for accessing protein-protein interaction databases.

Classes implementing this interface are:

  - :class:`BioGRID` for accessing `BioGRID <http://thebiogrid.org>`_
  - :class:`STRING` for accessing `CC`_ licensed `STRING <http://www.string-db.org/>`_
  - :class:`STRINGDetailed` for accessing `CC-NC-SA`_ licensed `STRING <http://www.string-db.org/>`_

.. _`CC`: http://creativecommons.org/licenses/by/3.0/

.. _`CC-NC-SA`: http://creativecommons.org/licenses/by-nc-sa/3.0/

The common interface
---------------------

.. autoclass:: PPIDatabase
   :members:
   :member-order: bysource

PPI databases
------------------

.. autoclass:: BioGRID
   :members:
   :inherited-members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: STRING
   :members:
   :inherited-members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: STRINGDetailed
   :members:
   :inherited-members:
   :member-order: bysource
   :show-inheritance:


