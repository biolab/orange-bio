import obiOMIM

diseases = obiOMIM.diseases()
genes = [obiOMIM.disease_genes(disease) for disease in diseases]

vertices = []
edges = []
for i in range(len(diseases)):
    vertices.append('%i "%s"\n' % (i + 1, diseases[i].name))
    for j in range(i + 1, len(diseases)):
        if set(genes[i]).intersection(genes[j]):
            edges.append('%i %i "%s"\n' %(i + 1, j + 1))

file = open("disease.net", "wb")
file.write("*Vertices %i 2\n" % len(vertices))
file.writelines(vertices)
file.write("*Edges\n")
file.writelines(edges)