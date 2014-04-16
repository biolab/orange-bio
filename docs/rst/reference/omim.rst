.. py:currentmodule:: Orange.bio.omim

.. index:: NCBI
.. index:: OMIM

**************************************************************
OMIM: Online Mendelian Inheritance in Man (:mod:`omim`)
**************************************************************

:obj:`omim` provides an interface to 
`Online Mendelian Inheritance in Man <http://www.ncbi.nlm.nih.gov/omim>`_
database. For now it only supports mapping genes to diseases.

.. autofunction:: diseases

.. autofunction:: genes

.. autofunction:: disease_genes

.. autofunction:: gene_diseases

The following example creates a network of connected diseases and
save it in `Pajek <http://vlado.fmf.uni-lj.si/pub/networks/pajek/>`_
.net format sets information file (:download:`omim_disease_network.py
<code/omim_disease_network.py>`).

.. literalinclude:: code/omim_disease_network.py

