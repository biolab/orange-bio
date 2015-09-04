============================================================
Dictyostelium discoideum databases (:mod:`dicty`)
============================================================

.. py:currentmodule:: orangecontrib.bio.dicty

The following example downloads experiments from the PIPA
database, specifically "RPKM + mapability expression (polyA) - R6"
results for all public experiments on Dictyostelium discoideum (dd)
at time point 16.

.. literalinclude:: code/pipax1.py

PIPAx database
--------------

.. autoclass:: orangecontrib.bio.dicty.PIPAx
   :members: __init__, genomes, get_data, mappings, result_types, results_list

DictyExpress database
---------------------

.. autoclass:: orangecontrib.bio.dicty.DictyExpress
   :members: __init__, objects, annotationTypes, annotations, get_data, search, annotationOptions

Auxillary functionality
-----------------------

.. autoclass:: orangecontrib.bio.dicty.CacheSQLite
    :members: __init__, list, contains, add, get, commit, clear

