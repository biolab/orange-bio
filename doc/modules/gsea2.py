import obiDicty
import obiGeneSets
import obiGsea
import orange

dbc = obiDicty.DatabaseConnection("http://asterix.fri.uni-lj.si/microarray/api/index.php?")
data = dbc.getData(sample='pkaC-', time="8")[0] #get first chip

print "First 10 examples"
for ex in data[:10]:
    print ex

genesets =  obiGeneSets.collections([":go:ddi", ":kegg:ddi"])
res = obiGsea.runGSEA(data, organism="ddi", minPart=0.05, geneSets=genesets, permutation="gene")

print "GSEA results"
print "%-40s %6s %6s %6s %7s" % ("LABEL", "NES", "P-VAL", "SIZE", "MATCHED") 
for name,resu in res.items()[:10]: 
    print "%-40s %6.3f %6.3f %6d %7d" % (name[:30], resu[1], resu[2], resu[4], resu[5]) 

