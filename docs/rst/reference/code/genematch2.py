import Orange.bio.kegg
import Orange.bio.gene

targets = Orange.bio.kegg.KEGGOrganism("9606").get_genes() #KEGG gene IDs

gmkegg = Orange.bio.gene.GMKEGG("9606")
gmgo = Orange.bio.gene.GMGO("9606")
gmkegggo = Orange.bio.gene.matcher([[gmkegg, gmgo]], direct=False) #joined matchers

gmkegg.set_targets(targets)
gmgo.set_targets(targets)
gmkegggo.set_targets(targets)

genes = [ "cct7", "pls1", "gdi1", "nfkb2", "a2a299" ]

print "%12s" % "gene", "%12s" % "KEGG", "%12s" % "GO", "%12s" % "KEGG+GO"
for gene in genes:
    print "%12s" % gene, "%12s" % gmkegg.umatch(gene), \
          "%12s" % gmgo.umatch(gene), \
          "%12s" % gmkegggo.umatch(gene)
