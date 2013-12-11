.. py:currentmodule:: Orange.bio.obiGene

.. index:: gene matching
.. index:: gene name matching
.. index:: matching
.. index:: NCBI
.. index:: gene info

********************************************************
Gene name matching (:mod:`obiGene`)
********************************************************

To use gene matchers
first set the target gene names with :obj:`~Matcher.set_targets` and then
match  with :obj:`~Matcher.match` or :obj:`~Matcher.umatch` functions. The
following example (:download:`genematch1.py <code/genematch1.py>`)
matches gene names to NCBI gene IDs:

.. literalinclude:: code/genematch1.py

Gene name matching
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


