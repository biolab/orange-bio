import Orange
from Orange.bio import dicty, geneset, gsea, gene

dbc = dicty.DictyExpress()
data = dbc.get_data(sample='pkaC-', time="8")

#select the first chip (the first attribute)
data = data.translate([data.domain.attributes[0]], True)

matcher = gene.matcher([[gene.GMKEGG("dicty"), gene.GMDicty()]])
genesets =  geneset.collections((("KEGG",), "dicty"))

res = gsea.direct(data, matcher=matcher, min_part=0.05, 
    gene_sets=genesets)

print "%-40s %6s %6s %6s %7s" % ("LABEL", "NES", "P-VAL", "SIZE", "MATCHED")
for name,resu in sorted(res.items()[:10], key=lambda x: x[1]["p"]): 
    print "%-40s %6.3f %6.3f %6d %7d" % (name.name[:35], resu["nes"],
        resu["p"], resu["size"], resu["matched_size"]) 

