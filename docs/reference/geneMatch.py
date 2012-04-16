import obiGene
import obiKEGG

targets = obiKEGG.KEGGOrganism("9606").get_genes() #human NCBI ID

gmkegg = obiGene.GMKEGG("9606")
gmgo = obiGene.GMGO("9606")
gmkegggo = obiGene.matcher([[gmkegg, gmgo]], direct=False)

gmkegg.set_targets(targets)
gmgo.set_targets(targets)
gmkegggo.set_targets(targets)

genes = [ "cct7", "pls1", "gdi1", "nfkb2", "dlg7" ]

print "%12s" % "gene", "%12s" % "KEGG", "%12s" % "GO", "%12s" % "KEGG+GO"
for gene in genes:
    print "%12s" % gene, "%12s" % gmkegg.umatch(gene), \
          "%12s" % gmgo.umatch(gene), \
          "%12s" % gmkegggo.umatch(gene)

