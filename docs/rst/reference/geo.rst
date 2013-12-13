.. py:currentmodule:: Orange.bio.obiGEO

.. index:: NCBI
.. index:: GEO
.. index:: Gene Expression Omnibus
.. index:: microarray data sets

**************************************************************
NCBI's Gene Expression Omnibus interface (:mod:`obiGEO`)
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

GDS classes
===========

.. autoclass:: GDSInfo
   :members:

An example that uses obj:`GDSInfo`::

    >>> from Orage import obiGEO
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
========

The following script prints out some information about a specific data
set. It does not download the data set, just uses the (local) GEO data
sets information file (:download:`geo_gds1.py <code/geo_gds1.py>`).

.. literalinclude:: code/geo_gds1.py
   :lines: 6-

The output of this script is::

    ID: GDS10
    Features: 39114
    Genes: 19883
    Organism: Mus musculus
    PubMed ID: 11827943
    Sample types:
      disease state (diabetic, diabetic-resistant, nondiabetic)
      strain (NOD, Idd3, Idd5, Idd3+Idd5, Idd9, B10.H2g7, B10.H2g7 Idd3)
      tissue (spleen, thymus)

    Description:
    Examination of spleen and thymus of type 1 diabetes nonobese diabetic
    (NOD) mouse, four NOD-derived diabetes-resistant congenic strains and
    two nondiabetic control strains.


GEO data sets provide a sort of mini ontology for sample labeling. Samples
belong to sample subsets, which in turn belong to specific types. Like
above GDS10, which has three sample types, of which the subsets for
the tissue type are spleen and thymus. For supervised data mining it
would be useful to find out which data sets provide enough samples for
each label. It is (semantically) convenient to perform classification
within sample subsets of the same type. We therefore need a script
that goes through the entire set of data sets and finds those, where
there are enough samples within each of the subsets for a specific
type. The following script does the work. The function ``valid``
determines which subset types (if any) satisfy our criteria. The
number of requested samples in the subset is by default set to ``n=40``
(:download:`geo_gds5.py <code/geo_gds5.py>`).

.. literalinclude:: code/geo_gds5.py
   :lines: 8-

The requested number of samples, ``n=40``, seems to be a quite
a stringent criteria met - at the time of writing of this documentation -
by 35 sample subsets. The output starts with::

    GDS1611
      genotype/variation: wild type/48, upf1 null mutant/48
    GDS3553
      cell type: macrophage/48, monocyte/48
    GDS3953
      protocol: training set/46, validation set/47
    GDS3704
      protocol: PUFA consumption/42, SFA consumption/42
    GDS3890
      agent: vehicle, control/46, TCE/48
    GDS1490
      other: non-neural/50, neural/100
    GDS3622
      genotype/variation: wild type/56, Nrf2 null/54
    GDS3715
      agent: untreated/55, insulin/55

Let us now pick data set GDS2960 and see if we can predict the disease
state. We will use logistic regression, and within 10-fold cross
validation measure AUC, the area under ROC. AUC is the probably for
correctly distinguishing between two classes if picking the sample from
target (e.g., the disease) and non-target class (e.g., control). From
(:download:`geo_gds6.py <code/geo_gds6.py>`)

.. literalinclude:: code/geo_gds6.py

The output of this script is::
    
    Samples: 101, Genes: 4068
    AUC = 0.960

The AUC for this data set is very high, indicating that using these gene
expression data it is almost trivial to separate the two classes.
