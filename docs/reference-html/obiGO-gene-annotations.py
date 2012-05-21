import obiGO

ontology = obiGO.Ontology()
annotations = obiGO.Annotations("sgd", ontology=ontology)

gene = annotations.aliasMapper["YIL075C"]
print gene, "(YIL075C) directly annotated to the folowing terms:"
for a in annotations.geneAnnotations[gene]:
    print ontology[a.GO_ID].name, "with evidence code", a.Evidence_code
    
# Get all genes annotated to the same terms as YIL075C
ids = set([a.GO_ID for a in annotations.geneAnnotations[gene]])
for GOID in ids:
	ants = annotations.GetAllAnnotations(GOID)
	genes = set([a.geneName for a in ants])
	print ", ".join(genes), "annotated to", GOID, ontology[a.GO_ID].name
	
	
