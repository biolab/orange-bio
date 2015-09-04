===================================
Array Express (:mod:`arrayexpress`)
===================================

.. py:currentmodule:: Orange.bio.arrayexpress

Access the `ArrayExpress`_ web services and database.

`ArrayExpress`_ is a database of gene expression experiments
that you can query and download.

Retrieve the object representing experiment with accession E-TABM-25

>>> from Orange.bio import arrayexpress
>>> experiment = ArrayExpressExperiment("E-TABM-25")
>>> print experiment.accession
E-TABM-25

>>> print experiment.name
Transcription profiling of aging in the primate brain

>>> print experiment.species
['Pan troglodytes']

>>> print experiment.files
[{'kind': ...

>>> # Retrieve the data matrix for experiment 'E-MEXP-2917'
>>> experiment = ArrayExpressExperiment("E-MEXP-2917")
>>> table = experiment.fgem_to_table()


Low level Array Express query using REST services:

>>> from Orange.bio import arrayexpress
>>> arrayexpress.query_experiments(accession='E-MEXP-31')
{u'experiments': ...

>>> arrayexpress.query_experiments(keywords='gliobastoma')
{u'experiments': ...

>>> arrayexpress.query_files(accession='E-MEXP-32', format="xml")
<xml.etree.ElementTree.ElementTree object ...

.. note:: Currently querying ArrayExpress files only works with the xml format.

.. _`ArrayExpress`: http://www.ebi.ac.uk/arrayexpress/


Interface
---------

.. autoclass:: ArrayExpressConnection
   :members:

.. autoclass:: ArrayDesign
   :members:

.. autoclass:: SampleDataRelationship
   :members:

.. autoclass:: InvestigationDesign
   :members:


Low-level querying with REST
-----------------------------

.. autofunction:: query_experiments

.. autofunction:: query_files


