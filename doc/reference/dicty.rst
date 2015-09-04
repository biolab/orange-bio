============================================================
Dictyostelium discoideum databases (:mod:`dicty`)
============================================================

.. py:currentmodule:: Orange.bio.dicty

The following example downloads experiments from the PIPA
database, specifically "RPKM + mapability expression (polyA) - R6"
results for all public experiments on Dictyostelium discoideum (dd)
at time point 16.

.. literalinclude:: code/pipax1.py

PIPAx database
--------------

.. autoclass:: Orange.bio.dicty.PIPAx
   :members: __init__, genomes, get_data, mappings, result_types, results_list

DictyExpress database
---------------------

.. autoclass:: Orange.bio.dicty.DictyExpress
   :members: __init__, objects, annotationTypes, annotations, get_data, search, annotationOptions

Auxillary functionality
-----------------------

.. autoclass:: Orange.bio.dicty.CacheSQLite
    :members: __init__, list, contains, add, get, commit, clear

