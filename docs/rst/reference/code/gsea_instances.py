import Orange
from Orange.bio import dicty, geneset, gsea, gene, geo

gds = geo.GDS("GDS10")
data = gds.getdata(transpose=True) 

matcher = gene.matcher([gene.GMKEGG("Homo sapiens")])
genesets = geneset.collections((("KEGG",), "Homo sapiens"))

#the number of permutations (n) should be much higher
res = gsea.run(data, gene_sets=genesets, matcher=matcher, 
    min_part=0.05, permutation="phenotype", n=10, 
    phen_desc=data.domain["tissue"], gene_desc=True) 

print
print "GSEA results (descriptor: tissue)"
print "%-40s %6s %6s %6s %7s" % ("LABEL", "NES", "FDR", "SIZE", "MATCHED") 
for gs, resu in sorted(res.items(), key=lambda x: x[1]["fdr"])[:10]: 
    print "%-40s %6.3f %6.3f %6d %7d" % (gs.name[:30],
        resu["nes"], resu["fdr"], resu["size"], resu["matched_size"]) 
