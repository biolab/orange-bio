"""
Print out some information on specific GEO's data set.
Does not download the data set.
"""

import orangecontrib.bio.geo
import textwrap

gdsinfo = orangecontrib.bio.geo.GDSInfo()
gds = gdsinfo["GDS10"]

print("ID:")
print(gds["dataset_id"])
print("Features: ")
print(gds["feature_count"])
print("Genes:")
print(gds["gene_count"])
print("Organism:")
print(gds["platform_organism"])
print("PubMed ID:")
print(gds["pubmed_id"])
print("Sample types:")
for sampletype in set([sinfo["type"] for sinfo in gds["subsets"]]):
    ss = [sinfo["description"] for sinfo in gds["subsets"] if sinfo["type"]==sampletype]
    print("  %s (%s)" % (sampletype, ", ".join(ss)))
print("")
print("Description:")
print("\n".join(textwrap.wrap(gds["description"], 70)))
