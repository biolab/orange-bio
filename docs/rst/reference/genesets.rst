.. py:currentmodule:: Orange.bio.obiGeneSets

.. index:: gene set
.. index:: gene sets

********************************************************
Gene sets (:mod:`obiGeneSets`)
********************************************************

:mod:`obiGeneSets` can load gene sets distributed with Orange (through ServerFiles) 
and also in `GMT file format <http://www.molmine.com/magma/fileformats.html>`_.

Load gene sets from KEGG and GO for mouse::

    obiGeneSets.collections((("KEGG",), "10090"), (("GO",), "10090"))

Open gene sets from ``specific.gmt`` file in the current working directory::

    obiGeneSets.collections("specific.gmt")

The above examples combined::

    obiGeneSets.collections((("KEGG",), "10090"), (("GO",), "10090"), "specific.gmt")


Loading gene sets
=================

.. autofunction:: list_all

.. autofunction:: collections


Supporting functionality
========================

.. autoclass:: GeneSets
   :members:

.. autoclass:: GeneSet
   :members:

.. autofunction:: register

