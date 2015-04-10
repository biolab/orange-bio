import orangecontrib.bio.kegg
import orangecontrib.bio.gene

targets = orangecontrib.bio.kegg.KEGGOrganism("9606").get_genes() #KEGG gene IDs

gmkegg = orangecontrib.bio.gene.GMKEGG("9606")
gmgo = orangecontrib.bio.gene.GMGO("9606")
gmkegggo = orangecontrib.bio.gene.matcher([[gmkegg, gmgo]], direct=False) #joined matchers

gmkegg.set_targets(targets)
gmgo.set_targets(targets)
gmkegggo.set_targets(targets)

genes = [ "cct7", "pls1", "gdi1", "nfkb2", "a2a299" ]

print("%12s %12s %12s %12s" % ( "gene", "KEGG", "GO", "KEGG+GO" ))
for gene in genes:
    print("%12s %12s %12s %12s" % \
        (gene, gmkegg.umatch(gene), gmgo.umatch(gene), gmkegggo.umatch(gene)))
