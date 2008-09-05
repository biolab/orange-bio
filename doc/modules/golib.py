import go
go.loadAnnotation("sgd")
go.loadGO()
geneNames=[ "YBR085W", ## synonim for AAC3
			"AAD10",
			"21S_rRNA_4" ##synonim for 21S_RRNA
			]
uniqueGeneNames = [go.geneMapper[name] for name in geneNames]
print uniqueGeneNames #should print ["AAC3", "AAD10", "21S_rRNA_4"]

#finding all terms with the keywords "catabolic" and "process" in their names
selectedTerms=filter(lambda term: "catabolic" in term.name and "process" in term.name, go.loadedGO.terms)
for term in selectedTerms:
	print term.id, term.name, term._def
	
def call(i): print i

res=go.GOTermFinder(uniqueGeneNames, progressCallback=call)
res=res.items()
res.sort(lambda (_1,(_2,a,_3)), (_4,(_5,b,_6)): cmp(a,b))
for GOId,(genes, p, n) in res:
	print GOId, p, genes

terms=go.GOTermFinder(uniqueGeneNames)
terms=go.filterByPValue(terms, 0.2)
go.drawEnrichmentGraph(terms, "enrichment_graph.png", width=600)