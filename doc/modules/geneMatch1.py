import obiGeneMatch
import obiKEGG

keggorg = obiKEGG.KEGGOrganism("mmu")
testgenes = keggorg.get_genes() 

targets = [ "Fndc4", "Itgb8", "Cdc34", "Olfr1403" ] 

gm = obiGeneMatch.GeneMatch(targets, organism="mmu")
match = gm.match(testgenes)
print "Matches", match

pnames = keggorg.list_pathways()

for keggname, origname in match:
	pwys = keggorg.get_pathways_by_genes([keggname])
	print origname, "is in", [ pnames[p] for p in pwys ] 
