import sys
import go
try:
    if sys.argv[1]=="remove":
        sys.exit(0)
except:
    pass

import obiKEGG, obiGO, obiGenomicsUpdate

pkgUpdate = obiGenomicsUpdate.PKGUpdate("go", obiGO.Update())

print "Updating GO ontology"
pkgUpdate.UpdateOntology()

for org, name in [("goa_human", "Homo sapiens"), ("sgd", "Yeast")]:
    print "Updating GO anotations for", name
    pkgUpdate.UpdateAnnotation(org)

pkgUpdate = obiGenomicsUpdate.PKGUpdate("kegg", obiKEGG.Update())

print "Updating KEGG taxonomy"
pkgUpdate.UpdateTaxonomy()

print "Updating KEGG orthology"
pkgUpdate.UpdateOrthology()

print "Updating KEGG reference pathways"
pkgUpdate.UpdateReference()

for org, name in [("hsa", "Homo sapiens"), ("sce", "Yeast")]:
    print "Updating KEGG pathways for", name
    pkg.UpdateOrganism(org)

