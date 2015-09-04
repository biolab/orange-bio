import orangecontrib.bio.gene
import orangecontrib.bio.taxonomy
import Orange

data = Orange.data.Table("brown-selected")

geneinfo = orangecontrib.bio.gene.NCBIGeneInfo('4932')

genes = [str(ex["gene"]) for ex in data]

for gene in genes:
    mappedgene = orangecontrib.bio.gene.homology.homolog(geneinfo(gene).symbol, \
        '4932', '9606')
    print(gene, mappedgene)
