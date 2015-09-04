import orangecontrib.bio.kegg
import orangecontrib.bio.gene

keggorg = orangecontrib.bio.kegg.KEGGOrganism("mmu")
kegg_genes = keggorg.get_genes() 

query = [ "Fndc4", "Itgb8", "Cdc34", "Olfr1403" ] 

gm = orangecontrib.bio.gene.GMKEGG("mmu") #use KEGG aliases for gene matching
gm.set_targets(kegg_genes) #set KEGG gene aliases as targets

for name in query:
    match = gm.umatch(name)
    if match:
        pwys = keggorg.get_pathways_by_genes([match])
        print(name + " is in")
        pathways = [ orangecontrib.bio.kegg.KEGGPathway(p).title for p in pwys ]
        if pathways:
            for a in pathways:
                print(' ' + a)
        else:
            print('  /')
