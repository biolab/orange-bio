import obiGene
import obiKEGG

keggorg = obiKEGG.KEGGOrganism("mmu")
kegg_genes = keggorg.get_genes() 

query = [ "Fndc4", "Itgb8", "Cdc34", "Olfr1403" ] 

gm = obiGene.GMKEGG("mmu") #use KEGG aliases for gene matching
gm.set_targets(kegg_genes) #set KEGG gene aliases as targets

pnames = keggorg.list_pathways()

for name in query:
    match = gm.umatch(name) # matched kegg alias or None
    if match:
    	pwys = keggorg.get_pathways_by_genes([match])
        print name, "is in", [ pnames[p] for p in pwys ] 
