import Orange.bio.dicty

pipa = Orange.bio.dicty.PIPAx()

results = pipa.results_list("R6")
dd16 = [ (i,d) for i,d in results.items() if \
        d["tp"] == '16' and d["species_id"] == "dd" ]

#group similar experiments with sorting
dd16 = sorted(dd16, key=lambda x: (x[1]["treatment"], x[1]["replicate"]))

data = pipa.get_data([i for i,d in dd16], exclude_constant_labels=True)

for at in data.domain.attributes:
    print at.name, \
        "treatment:", at.attributes["treatment"], \
        "replicate:", at.attributes["replicate"]

print

for a in data[:10]:
    print a
