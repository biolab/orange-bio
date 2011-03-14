import obiArrayExpress
from pprint import pprint

# test the gene_atlas_summary
summary = obiArrayExpress.get_atlas_summary(["Kalrn", "Ptprd", "Mbp", "Cyfip2"], "Mus musculus")
pprint(summary)

# test query_atlas_simple
results = obiArrayExpress.query_atlas_simple(genes=["Kalrn", "Ptprd", "Mbp", "Cyfip2"], organism="Mus musculus", regulation="up", condition="brain")
pprint(results)

# test Atlas Conditions
gene_cond1 = obiArrayExpress.AtlasConditionGeneProperty("Goterm", "Is", "translation")
gene_cond2 = obiArrayExpress.AtlasConditionGeneProperty("Disease", "Is", "cancer")
org_cond = obiArrayExpress.AtlasConditionOrganism("Homo sapiens")

conditions = obiArrayExpress.AtlasConditionList([gene_cond1, gene_cond2, org_cond])
results = obiArrayExpress.query_atlas(conditions)
pprint(results)

# test ArrayExpress experiments, files query

results = obiArrayExpress.query_experiments(accession="E-MEXP-31")
pprint(results)

results = obiArrayExpress.query_experiments(species="Homo sapines", expdesign="dose+response", ef="CellType")
pprint(results)

results = obiArrayExpress.query_experiments(species="Homo sapiens", gxa=True, assaycount=(1, 5), miamescore=(3, 5))
pprint(results)

results = obiArrayExpress.query_files(species="Mus musculus", gxa=True, keywords=["lung", "cancer"], miamescore=(3, 5), format="xml")
print results
