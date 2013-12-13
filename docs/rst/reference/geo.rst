.. py:currentmodule:: Orange.bio.obiGEO

.. index:: NCBI
.. index:: GEO
.. index:: Gene Expression Omnibus
.. index:: microarray data sets

**************************************************************
An interface to NCBI's Gene Expression Omnibus (:mod:`obiGEO`)
**************************************************************

:obj:`obiGEO` provides an interface to `NCBI
<http://www.ncbi.nlm.nih.gov/>`_'s `Gene Expression Omnibus
<http://www.ncbi.nlm.nih.gov/geo/>`_ repository. Currently, it only
supports `GEO DataSets <http://www.ncbi.nlm.nih.gov/sites/GDSbrowser>`_
information querying and retrieval.

The following illustrates how :obj:`GDS.getdata` is used to
construct a data set with genes in rows and samples in
columns. Notice that the annotation about each sample is retained
in ``.attributes``.

::

    >>> from Orange.bio import obiGEO
    >>> gds = obiGEO.GDS("GDS1676")
    >>> data = gds.getdata()
    >>> len(data)
    667
    >>> data[0]
    [?, ?, -0.803, 0.128, 0.110, -2.000, -1.000, -0.358], {"gene":'EXO1'}
    >>> data.domain.attributes[0]
    FloatVariable 'GSM63816'
    >>> data.domain.attributes[0].attributes
    Out[191]: {'dose': '20 U/ml IL-2', 'infection': 'acute ', 'time': '1 d'}

GDS
===

.. autoclass:: GDSInfo
   :members:

::

    >>> import obiGEO
    >>> info = obiGEO.GDSInfo()
    >>> info.keys()[:5]
    >>> ['GDS2526', 'GDS2524', 'GDS2525', 'GDS2522', 'GDS1618']
    >>> info['GDS2526']['title']
    'c-MYC depletion effect on carcinoma cell lines'
    >>> info['GDS2526']['platform_organism']
    'Homo sapiens'

.. autoclass:: GDS
   :members:


Examples
==================

Genes can have multiple aliases. When we combine data from different
sources, for example expression data with GO gene sets, we have to
match gene aliases representing the same genes. All implemented matching
methods are based on sets of gene aliases for one gene.

.. autoclass:: Matcher
   :members:

This modules provides the following gene matchers:

.. autoclass:: MatcherAliasesKEGG

.. autoclass:: MatcherAliasesGO

.. autoclass:: MatcherAliasesDictyBase

.. autoclass:: MatcherAliasesNCBI

.. autoclass:: MatcherAliasesEnsembl

.. autoclass:: MatcherDirect

Gene name matchers can applied in sequence (until the first match) or combined (overlapping sets of gene aliases of multiple gene matchers are combined) with the :obj:`matcher` function.

.. autofunction:: matcher

The following example tries to match input genes onto KEGG gene aliases (:download:`genematch2.py <code/genematch2.py>`).

.. literalinclude:: code/genematch2.py

Results show that GO aliases can not match onto KEGG gene IDs. For the last gene only joined GO and KEGG aliases produce a match::

        gene         KEGG           GO      KEGG+GO
        cct7    hsa:10574         None    hsa:10574
        pls1     hsa:5357         None     hsa:5357
        gdi1     hsa:2664         None     hsa:2664
       nfkb2     hsa:4791         None     hsa:4791
      a2a299         None         None     hsa:7052


The following example finds KEGG pathways with given genes (:download:`genematch_path.py <code/genematch_path.py>`).

.. literalinclude:: code/genematch_path.py

Output::

    Fndc4 is in
      /
    Itgb8 is in
      PI3K-Akt signaling pathway
      Focal adhesion
      ECM-receptor interaction
      Cell adhesion molecules (CAMs)
      Regulation of actin cytoskeleton
      Hypertrophic cardiomyopathy (HCM)
      Arrhythmogenic right ventricular cardiomyopathy (ARVC)
      Dilated cardiomyopathy
    Cdc34 is in
      Ubiquitin mediated proteolysis
      Herpes simplex infection
    Olfr1403 is in
      Olfactory transduction


