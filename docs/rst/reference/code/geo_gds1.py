"""
Print out some information on specific GEO's data set.
Does not download the data set.
"""

from Orange.bio import obiGEO
import textwrap

gdsinfo = obiGEO.GDSInfo()
gds = gdsinfo["GDS10"]

print "ID:", gds["dataset_id"]
print "Features:", gds["feature_count"]
print "Genes:", gds["gene_count"]
print "Organism:", gds["platform_organism"]
print "PubMed ID:", gds["pubmed_id"]
print "Sample types:"
for sampletype in set([sinfo["type"] for sinfo in gds["subsets"]]):
    ss = [sinfo["description"] for sinfo in gds["subsets"] if sinfo["type"]==sampletype]
    print "  %s (%s)" % (sampletype, ", ".join(ss))
print
print "Description:"
print "\n".join(textwrap.wrap(gds["description"], 70))
