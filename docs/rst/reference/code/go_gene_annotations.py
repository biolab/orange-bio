from Orange.bio import go

ontology = go.Ontology()

# Print names and definitions of all terms with "apoptosis" in the name
apoptosis = [term for term in ontology.terms.values()
             if "apoptosis" in term.name.lower()]
for term in apoptosis:
    print(term.name + term.id)
    print(term.def_)

# Load annotations for yeast.
annotations = go.Annotations("sgd", ontology=ontology)

res = annotations.get_enriched_terms(["YGR270W", "YIL075C", "YDL007W"])

gene = annotations.alias_mapper["YIL075C"]
print(gene + " (YIL075C) directly annotated to the following terms:")
for a in annotations.gene_annotations[gene]:
    print(ontology[a.GO_ID].name + " with evidence code " + a.Evidence_Code)

# Get all genes annotated to the same terms as YIL075C
ids = set([a.GO_ID for a in annotations.gene_annotations[gene]])
for termid in ids:
    ants = annotations.get_all_annotations(termid)
    genes = set([a.DB_Object_Symbol for a in ants])
    print(", ".join(genes) +" annotated to " + termid + " " + ontology[termid].name)
