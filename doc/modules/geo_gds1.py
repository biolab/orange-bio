import obiGEO

gds_info = obiGEO.GDSInfo()

gds = gds_info["GDS10"]

print "ID:", gds["dataset_id"]
print "Features:", gds["feature_count"]
print "Genes:", gds["gene_count"]
print "Organism:", gds["platform_organism"]
print "PubMed ID:", gds["pubmed_id"]
print "Sample types:", ", ".join(set([sampleinfo["type"] for sampleinfo in gds["subsets"]]))

print
print "Description"
print gds["description"]
