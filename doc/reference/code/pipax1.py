import orangecontrib.bio.dicty

pipa = orangecontrib.bio.dicty.PIPAx()

results = pipa.results_list("R6")
dd16 = [ (i,d) for i,d in results.items() if \
        d["tp"] == '16' and d["species_id"] == "dd" ]

#group similar experiments with sorting
dd16 = sorted(dd16, key=lambda x: (x[1]["treatment"], x[1]["replicate"]))

data = pipa.get_data([i for i,d in dd16], exclude_constant_labels=True, \
    allowed_labels=["id", "treatment", "replicate"])

def print_data(data):
    for at in data.domain.attributes:
        print("%s treatment: %s replicate: %s" % \
            (at.name, at.attributes["treatment"], at.attributes["replicate"]))
    print("")
    for a in data[:10]:
        print(a)

print_data(data)

print("")
datar = orangecontrib.bio.dicty.join_replicates(data)

print_data(datar)
