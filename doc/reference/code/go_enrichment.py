from orangecontrib.bio import go

ontology = go.Ontology()
annotations = go.Annotations("sgd", ontology=ontology)

res = annotations.get_enriched_terms(["YGR270W", "YIL075C", "YDL007W"])
print("Enriched terms:")
for go_id, (genes, p_value, ref) in res.items():
    if p_value < 0.05:
        print(ontology[go_id].name + " with p-value: %.4f " % p_value + ", ".join(genes))

# And again for slims
ontology.set_slims_subset("goslim_yeast")

res = annotations.get_enriched_terms(["YGR270W", "YIL075C", "YDL007W"],
                                     slims_only=True)
print("Enriched slim terms:")
for go_id, (genes, p_value, _) in res.items():
    if p_value < 0.05:
        print(ontology[go_id].name + " with p-value: %.4f " % p_value + ", ".join(genes))
