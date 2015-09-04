.. py:currentmodule:: orangecontrib.bio.geneset
.. py:module:: orangecontrib.bio.geneset

.. index:: gene set
.. index:: gene sets

********************************************************
Gene sets (:mod:`geneset`)
********************************************************

This module can load either gene sets distributed with Orange
or custom gene sets in the `GMT file format <http://www.molmine.com/magma/fileformats.html>`_.

The available gene set collection can be listed with :obj:`list_all`.

:obj:`collections` loads gene sets. Gene sets provided with Orange are
organized hierarchically. Although the ``GO`` hierarchy includes subsets,
all of them can be loaded with (the organism here is a mouse)::

    orangecontrib.bio.geneset.collections((("GO",), "10090"))

To open multiple gene set collections at once, for example, KEGG and GO, try::

    orangecontrib.bio.geneset.collections((("KEGG",), "10090"), (("GO",), "10090"))

You could also open a file with gene sets. The following line would open
``specific.gmt`` from the current working directory::

    orangecontrib.bio.geneset.collections("specific.gmt")

The above examples combined::

    orangecontrib.bio.geneset.collections((("KEGG",), "10090"), (("GO",), "10090"), "specific.gmt")

Furthermore, all gene sets for a specific organism can be opened with an empty
hierarchy::

    orangecontrib.bio.geneset.collections((tuple(), "10090"))


Loading gene sets
=================

.. autofunction:: orangecontrib.bio.geneset.list_all

.. autofunction:: orangecontrib.bio.geneset.collections


Supporting functionality
========================

.. autoclass:: orangecontrib.bio.geneset.GeneSets
   :members:
   :show-inheritance:

.. autoclass:: orangecontrib.bio.geneset.GeneSet
   :members:

.. autofunction:: orangecontrib.bio.geneset.register

