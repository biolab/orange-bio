import orange
import obiHomoloGene, obiTaxonomy, obiGene

data = orange.ExampleTable("../../../../doc/datasets/brown-selected")

geneinfo = obiGene.NCBIGeneInfo('4932')

genes = [str(ex["gene"]) for ex in data]

for gene in genes:
    mappedgene = obiHomoloGene.homolog(geneinfo(gene).symbol, '4932', '9606')
    print gene, mappedgene