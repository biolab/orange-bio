import obiGO

ontology = obiGO.Ontology.Load()
annotations = obiGO.Annotations.Load("sgd", ontology=ontology)

ontology.SetSlimsSubset("goslim_yeast")
terms = annotations.GetAnnotatedTerms(["YGR270W", "YIL075C", "YDL007W"], directAnnotationOnly=True)
slims = set()
for term in terms:
    print term
    slims.update(ontology.GetSlimTerms(term))

print "Genes: YGR270W, YIL075C and YDL007W map to the folowing slims terms:"
for term in slims:
    print term, ontology.terms[term].name