"""
Getting genesets from KEGG and GO.

Maintainer: Marko Toplak
"""

import obiKEGG, obiGeneMatch, orange

import go

def allparentsi(goid, parentsd):
    parents = go.loadedGO.termDict[goid].parents
    parentsr =  parents + reduce(lambda x,y: x+y, [ allparents(pid, parentsd) for pid in parents ], [])
    return parentsr

def allparents(goid, parentsd=None):
    if parentsd != None and goid in parentsd:
        return parentsd[goid]
    else:
        parentsr = list(set(allparentsi(goid, parentsd)))
        if parentsd != None:
            parentsd[goid] = parentsr
        return parentsr

parentsd = {}

def _geneToTerms(goorg):
    go.loadGO()
    go.loadAnnotation(goorg)

    gt = {}
    for gene,data in go.loadedAnnotation.geneAnnotations.items():
        for ant in data:
            els = gt.get(gene,set())
            els.add(ant.GOId)
            gt[gene] = els

    directAnnotationOnly = False
    if not directAnnotationOnly:
        for gene, data in gt.items():
            addp = set()
            for el in data:
                addp.update(set(allparents(el, parentsd)))
            data.update(addp)
            gt[gene] = data

    for gene, data in gt.items():
        gt[gene] = sorted(data)    
 
    return gt

def goGeneSets(goorg):
    """
    Returns gene sets from GO. Look at the annotation
    of the organism provided.
    Returns a list of [ (GOid, genes) ]
    """
    go.loadGO()

    gt = _geneToTerms(goorg)

    td = go.loadedGO.termDict

    map = {} #map from set id to to it's genes

    for gene,sets in gt.items():
        for s in sets:
            cs = map.get(s, set())
            cs.add(gene)
            map[s] = cs

    for a,b in map.items():
        gt[a] = sorted(b)
 
    return map

def keggGeneSets(keggorg):
    """
    Returns pathways from KEGG for provided organism.
    Returns a list of [ (name, genes) ]
    """
    kegg = obiKEGG.KEGGOrganism(keggorg)

    pways = kegg.list_pathways()

    dicp = {}
    for id,desc in pways.items():
        dicp[desc] = kegg.get_genes_by_pathway(id)

    return dicp

def addSource(dic, addition):
    return dict( \
        [ (addition + name,genes) for name,genes in dic.items() ]\
        )

if __name__ == "__main__":
    print keggGeneSets("sce").items()[:10]
    print goGeneSets("sgd").items()[:10]

