from Orange.bio import obiGene

#matching targets are NCBI gene IDs
targets = obiGene.NCBIGeneInfo("Homo sapiens").keys()

gm = obiGene.GMNCBI("9606")
gm.set_targets(targets)

for gene in [ "cct7", "pls1", "gdi1", "nfkb2", "dlg7" ]:
    print 'Gene', gene, 'is NCBI gene', gm.umatch(gene)
