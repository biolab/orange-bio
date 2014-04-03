import Orange.bio.omim

diseases = Orange.bio.omim.diseases()
genes = [Orange.bio.omim.disease_genes(disease) for disease in diseases]

vertices = []
edges = []
for i in range(len(diseases)):
    vertices.append('%i "%s"\n' % (i + 1, diseases[i].name))
    for j in range(i + 1, len(diseases)):
        intersection = set(genes[i]).intersection(genes[j])
        if intersection:
            edges.append('%i %i %i l "%s"\n' %(i + 1, j + 1, len(intersection), ",".join(sorted(intersection))))

file = open("disease.net", "wb")
file.write("*Vertices %i 2\n" % len(vertices))
file.writelines(vertices)
file.write("*Edges\n")
file.writelines(edges)
