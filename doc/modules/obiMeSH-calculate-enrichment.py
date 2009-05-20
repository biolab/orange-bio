import obiMeSH
import orange

# load datasets
reference = orange.ExampleTable('obiMeSH-reference-dataset.tab')
cluster = orange.ExampleTable('obiMeSH-cluster-dataset.tab')

# find and print enriched MeSH terms with p-value < 0.1
d = obiMeSH.obiMeSH()
enrichment = d.findEnrichedTerms(reference, cluster, pThreshold=0.1)
d.printMeSH(enrichment)