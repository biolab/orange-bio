"""
Getting genesets from KEGG and GO.

Maintainer: Marko Toplak
"""

import obiKEGG, obiGeneMatch, orange
import os
import go

def nth(l,n):
    return [ a[n] for a in l]

collectionsPath = None

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
    Returns a ditionary  of (GOid, genes)
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

    nmap = {}
    for a,b in map.items():
        nmap[a] = sorted(b)
 
    return nmap

def keggGeneSets(keggorg):
    """
    Returns pathways from KEGG for provided organism.
    Returns a dictionary if (name, genes)
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

def keggToGo(org):
    """
    Converts KEGG organism to GO organism.
    """
    map = dict([ (b,a) for a,b in obiGeneMatch.GeneMatchMk2.dbOrgMap.items()])
    return map.get(org, None)

"""
CUSTOM GENESETS
"""

def strip(s):
    #return s
    return s.rstrip().lstrip()

def possiblyReadFile(s):
    """
    If s is not a string, then it is probably a file.
    Read it's contents.
    """
    if isinstance(s, basestring):
        return s
    else:
        return s.read()

def handleNELines(s, fn):
    """
    Run function on nonempty lines of a string.
    Return a list of results for each line.
    """
    lines = s.split("\n")
    lines = [ strip(l) for l in lines ]
    lines = filter(lambda x: x != "", lines)
    return [ fn(l) for l in lines ]

def genesetsLoadGMT(s):
    """
    Eech line consists of tab separated elements. First is
    the geneset name, next is it's description. 
    
    For now the description is skipped.
    """

    s = possiblyReadFile(s)

    def hline(s):
        tabs = [ strip(tab) for tab in s.split("\t") ]
        return tabs[0], tabs[2:]

    return dict(handleNELines(s, hline))

def collectionsPathname():
    if not collectionsPath:
        import orngEnviron
        pth = os.path.join(orngEnviron.directoryNames["bufferDir"], "gsea", "genesets")
        try:
            os.makedirs(pth)
        except:
            pass

        return pth
    else:
        return collectionsPath


def prettyfygo(god):
    import go
    go.loadGO()

    def translatens(a):
        if (a == "molecular_function"): return "MF"
        elif (a == "biological_process"): return "BP"
        elif (a == "cellular_component"): return "CC"
        else: 
            print a
            ffkkfds()

    ndic = {}
    for a,b in god.items():
        ndic["[GO " + translatens(go.loadedGO.termDict[a].namespace) + "] " + go.loadedGO.termDict[a].name] = b
    return ndic


def createCollection(lnf):
    """
    Input - list of tuples of geneset collection name and GMT
    file.
    """

    def addSource(dic, addition):
        return dict( \
            [ (addition + name,genes) for name,genes in dic.items() ] )

    gen1 = {}

    for n,fn in lnf:

        if fn.lower()[-4:] == ".gmt": #format from webpage
            gen2 = genesetsLoadGMT(open(fn,"rt"))
            gen1.update(addSource(gen2, "[%s] " % n))

        elif fn.lower()[-4:] == ".pck": #pickled dictionary
            import pickle
            f = open(fn,'rb')
            gen2 = pickle.load(f)
            gen1.update(addSource(gen2, "[%s] " % n))

        elif n == fn and n[0] == ":":
            _,name,org = n.split(":")

            if name == "kegg":
                gen2 = keggGeneSets(org)
                gen1.update(addSource(gen2, "[%s] " % "KEGG"))

            elif name == "go":
                gen2 = goGeneSets(keggToGo(org))
                gen2 = prettyfygo(gen2)
                gen1.update(addSource(gen2, ""))

            else:
                raise Exception("Wrong special name (%s)" % (name))            
            
        else:
            raise Exception("Can not recognize geneset (%s,%s)" % (n,fn))

    return gen1

def getCollectionFiles(path=collectionsPathname()):

    def loadInfo(path):
        #TODO load info file
        return {}
        
    info = loadInfo(path)

    def okFile(fn):
        if fn.lower()[-4:] == ".gmt":
            return True
        elif fn.lower()[-4:] == ".pck":
            return True
        return False

    files = sorted(filter(okFile, os.listdir(path)))

    #print files, os.listdir(path)

    out = []

    for file in files:
        fn = os.path.join(path, file)
        name = info.get(file, file)
        if name == file:
            #remove suffix if no info
            if name.lower()[-4:] == ".gmt" or name.lower()[-4:] == ".pck":
                name = name[:-4]

        out.append( (name, fn) )

    return out

def collections(l=[], default=True, path=collectionsPathname()):
    """
    Input is a list of collections.
    Default - if default collections are included.
    Input is a list of names. If names match to any names in path,
    they are taken. If not, file with that name is regarded as
    a filename of gene set colections
    """
    collections = getCollectionFiles(path) #default collections

    coln = nth(collections, 0)
    colff = nth(collections, 1)
    colf =  [ os.path.split(f)[1] for f in colff ]

    check = [ coln, colff, colf ]

    choosen = set()
    if default:
        choosen = choosen | set(collections)

    for col in l:
        added = False
        if not added:
            try: # if integer it can be the index
                choosen = choosen | set( [ collections[int(col)] ])
                added = True
            except:
                pass
        if not added:
            for cl in check:
                if col in cl:
                    choosen = choosen | set( [ collections[cl.index(col)] ])
                    added = True
                    break
        if not added:
            choosen = choosen | set( [ (col, col) ] )

    #pair in choosen are (name, location)

    return createCollection(list(choosen))

"""
End genesets
"""


if __name__ == "__main__":
    #print keggGeneSets("sce").items()[:10]
    #print goGeneSets("sgd").items()[:10]

    print keggToGo("mmu")
    print keggToGo2("mmu")



