import obiArrayExpress

# test the gene_atlas_summary
from pprint import pprint
summary = obiArrayExpress.get_atlas_summary(["Kalrn", "Ptprd", "Mbp", "Cyfip2"], "Mus musculus")
pprint(summary)

# test Atlas Conditions

gene_cond1 = obiArrayExpress.AtlasConditionGeneProperty("Goterm", "Is", "translation")
gene_cond2 = obiArrayExpress.AtlasConditionGeneProperty("Disease", "Is", "cancer")
org_cond = obiArrayExpress.AtlasConditionOrganism("Homo sapiens")

conditions = obiArrayExpress.AtlasConditionList([gene_cond1, gene_cond2, org_cond])
results = obiArrayExpress.query_atlas(conditions, format="json")
import json
results = json.load(results)
pprint(results)

# test ArrayExpress experiments, files query

#from xml.dom.minidom import parse
from json import load as parse
results = parse(obiArrayExpress.query_experiments(accession="E-MEXP-31", format="json"))
pprint(results)

results = parse(obiArrayExpress.query_experiments(species="Homo sapines", expdesign="dose+response", ef="CellType", format="json"))
pprint(results)



