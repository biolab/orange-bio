============================================================
Gene Set Enrichment Analysis (GSEA, :mod:`gsea`)
============================================================

.. py:currentmodule:: Orange.bio.gsea


Gene Set Enrichment Analysis (GSEA) [GSEA]_ aims to identify enriched gene sets
given gene expression data for multiple samples with their phenotypes.

.. autofunction:: Orange.bio.gsea.run

.. autofunction:: Orange.bio.gsea.direct

Examples: gene expression data
------------------------------

The following examples use a gene expression data set from the GEO database.
We show the same analysis on two formats of data.

With samples as instances (in rows):

.. literalinclude:: code/gsea_instances.py

With samples as features (in columns):

.. literalinclude:: code/gsea_genes.py

Both scripts output::

    GSEA results (descriptor: tissue)
    LABEL                                       NES    FDR   SIZE MATCHED
    Porphyrin and chlorophyll meta           -1.817  0.000     43      23
    Staphylococcus aureus infectio           -1.998  0.000     59      28
    Non-homologous end-joining                1.812  0.000     13      12
    Fanconi anemia pathway                    1.911  0.000     53      27
    Cell cycle                                1.777  0.000    124     106
    Glycine, serine and threonine            -2.068  0.000     39      29
    HIF-1 signaling pathway                  -1.746  0.000    106      90
    Ether lipid metabolism                   -1.788  0.000     42      27
    Fc epsilon RI signaling pathwa           -1.743  0.000     70      53
    B cell receptor signaling path           -1.782  0.000     72      62

Example: our own gene sets
--------------------------

We present a simple example on iris data set. Because data set
is not a gene expression data set, we had to specify our own
sets of features that belong together.

.. literalinclude:: code/gsea1.py

The output::

    LABEL     NES  P-VAL GENES
    sepal   1.087  0.630 ['sepal width', 'sepal length']
    petal  -1.117  0.771 ['petal width', 'petal length']


Example: directly passing correlation data
------------------------------------------

GSEA can also directly use correlation data between individual genes
and a phenotype. If (1) input data with only one example (attribute names 
are gene names) or (2) there is only one continuous feature in the given 
data set (gene names are in the first :obj:`Orange.feature.String`.
The following outputs ten pathways with smallest p-values.

.. literalinclude:: code/gsea2.py

The output::

    LABEL                                       NES  P-VAL   SIZE MATCHED
    Biosynthesis of amino acids               1.407  0.056     58      40
    beta-Alanine metabolism                   1.165  0.232     13      10
    Taurine and hypotaurine metabolism        1.160  0.413      4       3
    Porphyrin and chlorophyll metabolis      -0.990  0.517     14       5
    Valine, leucine and isoleucine degr       0.897  0.585     29      21
    Ether lipid metabolism                    0.713  0.857     10       6
    Biosynthesis of unsaturated fatty a       0.659  0.922     10       6
    Protein processing in endoplasmic r       0.647  0.941     71      40
    RNA polymerase                            0.550  0.943     24       7
    Glycosylphosphatidylinositol(GPI)-a      -0.540  0.946     19       4

.. [GSEA] Subramanian, Aravind et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 2005.
