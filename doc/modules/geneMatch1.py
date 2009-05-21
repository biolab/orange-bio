import obiGene
import obiKEGG

keggorg = obiKEGG.KEGGOrganism("mmu")
kegg_genes = keggorg.get_genes() 

targets = [ "Fndc4", "Itgb8", "Cdc34", "Olfr1403" ] 

gm = obiGene.GMKEGG("mmu") #use KEGG aliases for gene matching
gm.set_targets(targets)

pnames = keggorg.list_pathways()

for keggname in kegg_genes:
    match = gm.umatch(keggname)
    if match:
        pwys = keggorg.get_pathways_by_genes([keggname])
        print match, "is in", [ pnames[p] for p in pwys ] 
