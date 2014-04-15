=========================
Bio Mart (:mod:`biomart`)
=========================

.. py:currentmodule:: Orange.bio.biomart

Access BioMart MartService.

>>> from Orange.bio.biomart import *
>>> connection = BioMartConnection(
...     "http://www.biomart.org/biomart/martservice")
...
>>> reg = BioMartRegistry(connection)
>>> for mart in reg.marts():
...    print mart.name
...
ensembl...
>>> dataset = BioMartDataset(
...     mart="ensembl", internalName="hsapiens_gene_ensembl",
...     virtualSchema="default", connection=connection)
...
>>> for attr in dataset.attributes()[:10]:
...    print attr.name
...
Ensembl Gene ID...
>>> data = dataset.get_data(
...    attributes=["ensembl_gene_id", "ensembl_peptide_id"],
...    filters=[("chromosome_name", "1")])
...
>>> query = BioMartQuery(reg.connection, virtualSchema="default")
>>> query.set_dataset("hsapiens_gene_ensembl")
>>> query.add_attribute("ensembl_gene_id")
>>> query.add_attribute("ensembl_peptide_id")
>>> query.add_filter("chromosome_name", "1")
>>> count = query.get_count()

Interface
---------

.. autoclass:: BioMartConnection
   :members:

.. autoclass:: BioMartRegistry
   :members:

.. autoclass:: BioMartQuery
   :members:

.. autoclass:: BioMartDataset
   :members:

.. autoclass:: BioMartVirtualSchema
   :members:

.. autoclass:: BioMartDatabase
   :members:

.. autoclass:: Attribute
   :members:

.. autoclass:: Filter
   :members:

.. autoclass:: DatasetConfig
   :members:
