import Orange
from Orange.bio import obiGsea, obiGene

data = Orange.data.Table("iris")

gen1 = dict([
    ("sepal",["sepal length", "sepal width"]), 
    ("petal",["petal length", "petal width", "petal color"])
    ])

res = obiGsea.runGSEA(data, matcher=obiGene.matcher([]), minSize=2, geneSets=gen1)
print "%5s  %6s %6s %s" % ("LABEL", "NES", "P-VAL", "GENES")
for gs,resu in res.items():
    print "%5s  %6.3f %6.3f %s" % (gs.id, resu["nes"], resu["p"], str(resu["genes"]))
