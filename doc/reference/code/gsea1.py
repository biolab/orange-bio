import Orange
import orangecontrib.bio.gsea
import orangecontrib.bio.gene
import orangecontrib.bio.geneset

data = Orange.data.Table("iris")

gen1 = orangecontrib.bio.geneset.GeneSets(dict([
    ("sepal", ["sepal length", "sepal width"]), 
    ("petal", ["petal length", "petal width", "petal color"])
    ]))

res = orangecontrib.bio.gsea.run(data, gene_sets=gen1, matcher=orangecontrib.bio.gene.matcher([]), min_size=2)
print "%5s  %6s %6s %s" % ("LABEL", "NES", "P-VAL", "GENES")
for gs,resu in res.items():
    print "%5s  %6.3f %6.3f %s" % (gs.id, resu["nes"], resu["p"], str(resu["genes"]))
