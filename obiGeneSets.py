"""
Getting genesets from KEGG and GO.

Maintainer: Marko Toplak
"""

"""
TODO - add addtitional description to GENE SETS 
name,
desciption,
link
"""

import obiKEGG, orange
import os
import obiGO

def nth(l,n):
    return [ a[n] for a in l]

collectionsPath = None

class GeneSet(object):

    def __init__(self, pair=None, genes=None, name=None, id=None, \
        description=None, link=None):
        """
        pair can be (id, listofgenes) - it is used before anything else.
        """
        if pair:
            self.id, self.genes = pair[0], set(pair[1])
            self.name = self.id

        if genes == None:
            genes = []

        self.genes = set(genes)
        self.name = name
        self.id = id
        self.description = description
        self.link = link

    def odict(self):
        """
        Returns a pair (id, listofgenes), like in old format.
        """
        return self.id, self.genes

class GeneSetIDException(Exception):
    pass

class GeneSets(object):
    
    def init(self, odict=None, name=None, lgs=None):
        """
        odict are genesets in old dict format.
        lgs are 
        """
        self.name = name
        self.idict = {}
        if odict != None:
            self.idict = dict((i,GeneSet(pair=(i,g))) for i,g in odict.items())
        elif lgs != None:
            self.idict = dict((g.id,g) for g in lgs)

    def keys(self):
        return self.idict.keys()

    def get(self, a):
        return self.idict[a]

    def odict(self):
        """ Return gene sets in old dictionary format. """
        return dict(gs.odict() for gs in self.idict.values())

    def add(self, gs):
        if gs.id in idict:
            raise GeneSetIDException

def geneSetUnion():
    pass

def goGeneSets(goorg):
    """
    Returns gene sets from GO. Look at the annotation
    of the organism provided.
    Returns a ditionary  of (GOid, genes)
    """

    ontology = obiGO.Ontology.Load()
    annotations = obiGO.Annotations.Load(goorg, ontology=ontology)

    terms = ontology.terms.keys()

    map = {} #map from set id to to it's genes

    for term in terms:
        genes = annotations.GetAllGenes(term)
        if len(genes):
            map[term] = genes

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

    pways = kegg.pathways()

    dicp = {}
    for id in pways:
        pway = obiKEGG.KEGGPathway(id)
        dicp[pway.title] = kegg.get_genes_by_pathway(id)

    return dicp

def addSource(dic, addition):
    return dict( \
        [ (addition + name,genes) for name,genes in dic.items() ]\
        )

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

def loadGMT(s):
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

def loadPickled(s):
    import cPickle
    s = possiblyReadFile(s)
    gen2 = cPickle.loads(s)
    return gen2

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
    ontology = obiGO.Ontology.Load()
 
    def translatens(a):
        if (a == "molecular_function"): return "MF"
        elif (a == "biological_process"): return "BP"
        elif (a == "cellular_component"): return "CC"
        else: 
            print a
            ffkkfds()

    ndic = {}
    for a,b in god.items():
        ndic["[GO " + translatens(ontology.terms[a].namespace) + "] " + a + " " + ontology.terms[a].name] = b
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
            gen2 = loadGMT(open(fn,"rt"))
            gen1.update(addSource(gen2, "[%s] " % n))

        elif fn.lower()[-4:] == ".pck": #pickled dictionary
            gen2 = loadPickled(open(fn,"rb"))
            gen1.update(addSource(gen2, "[%s] " % n))

        elif n == fn and n[0] == ":":
            _,name,org = n.split(":")

            if name == "kegg":
                gen2 = keggGeneSets(org)
                gen1.update(addSource(gen2, "[%s] " % "KEGG"))

            elif name == "go":
                gen2 = goGeneSets(org)
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

def collections(l=[], default=False, path=collectionsPathname()):
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
    print keggGeneSets("sce").items()[:10]
    #col = collections([":go:sce"])
    #print col.items()[:10]
    
