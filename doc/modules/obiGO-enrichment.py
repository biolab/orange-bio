import obiGO

ontology = obiGO.Ontology.Load()
annotations = obiGO.Annotations.Load("sgd", ontology=ontology)

res = annotations.GetEnrichedTerms(["YGR270W", "YIL075C", "YDL007W"])
print "Enriched terms:"
for GOId, (genes, p_value, ref) in res.items():
    if p_value < 0.05:
        print ontology.terms[GOId].name, "with p-value: %.4f" %p_value, ", ".join(genes)

# And again for slims        
ontology.SetSlimsSubset("goslim_yeast")

res = annotations.GetEnrichedTerms(["YGR270W", "YIL075C", "YDL007W"], slimsOnly=True)
print "Enriched slim terms:"
for GOId, (genes, p_value, _) in res.items():
    if p_value < 0.05:
        print ontology.terms[GOId].name, "with p-value: %.4f" %p_value, ", ".join(genes)
        
# Print names and definitions of all terms with "apoptosis" in the name
##for term in [term for term in ontology.terms.values() if "apoptosis" in term.name.lower()]:
##	print term.name, term.id
##	print term.def_ if hasattr(term, "def_") else ""