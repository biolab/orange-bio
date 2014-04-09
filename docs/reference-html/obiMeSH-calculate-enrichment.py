from Orange.bio import obiMeSH
import Orange

# load datasets
reference = Orange.data.Table('obiMeSH-reference-dataset.tab')
cluster = Orange.data.Table('obiMeSH-cluster-dataset.tab')

# find and print enriched MeSH terms with p-value < 0.1
d = obiMeSH.obiMeSH()
enrichment = d.findEnrichedTerms(reference, cluster, pThreshold=0.1)
d.printMeSH(enrichment)
