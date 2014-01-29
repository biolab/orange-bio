import Orange
from Orange.bio import obiDicty, obiGeneSets, obiGsea,obiGene

dbc = obiDicty.DatabaseConnection()
data = dbc.get_single_data(sample='pkaC-', time="8")

#select the first chip (the first attribute)
data = data.translate([data.domain.attributes[0]], True)

matcher = obiGene.matcher([[obiGene.GMKEGG("dicty"), obiGene.GMDicty()]])
genesets =  obiGeneSets.collections((("KEGG",), "dicty"))

res = obiGsea.runGSEA(data, matcher=matcher, minPart=0.05, 
    geneSets=genesets, permutation="gene")

print "%-40s %6s %6s %6s %7s" % ("LABEL", "NES", "P-VAL", "SIZE", "MATCHED")
for name,resu in sorted(res.items()[:10], key=lambda x: x[1]["p"]): 
    print "%-40s %6.3f %6.3f %6d %7d" % (name.name[:35], resu["nes"],
        resu["p"], resu["size"], resu["matched_size"]) 

