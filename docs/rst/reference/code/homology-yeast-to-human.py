import Orange.bio.gene.homology
import Orange.bio.taxonomy

data = Orange.data.Table("brown-selected")

geneinfo = Orange.bio.gene.NCBIGeneInfo('4932')

genes = [str(ex["gene"]) for ex in data]

for gene in genes:
    mappedgene = Orange.bio.gene.homology.homolog(geneinfo(gene).symbol, \
        '4932', '9606')
    print gene, mappedgene
