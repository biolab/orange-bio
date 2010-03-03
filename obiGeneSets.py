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
        description=None, link=None, organism=None, hierarchy=None):
        """
        pair can be (id, listofgenes) - it is used before anything else.
        """
        if pair:
            self.id, self.genes = pair[0], set(pair[1])
            self.name = self.id

        if genes == None:
            genes = []

        self.hierarchy = hierarchy
        self.genes = set(genes)
        self.name = name
        self.id = id
        self.description = description
        self.link = link
        self.organism = organism

    def size(self):
        return len(self.genes)

    def to_odict(self, source=True, name=True):
        """
        Returns a pair (id, listofgenes), like in old format.
        """
        oname = self.id
        if source and self.hierarchy:
            oname = "[ " + ", ".join(self.hierarchy) + " ] " + oname
        if name and self.name:
            oname = oname + " " + self.name

        return oname, self.genes

    def __repr__(self):
        return "GeneSet(" + ", ".join( [ 
            "id=" + str(self.id),
            "genes=" + str(self.genes),
            "name=" + str(self.name),
            "link=" + str(self.link),
            "hierarchy=" + str(self.hierarchy)
        ]) + ")"

class GeneSetIDException(Exception):
    pass

class GeneSets(object):
    
    def __init__(self, odict=None, name=None, gs=None):
        """
        odict are genesets in old dict format.
        gs are genesets in new format
        """
        self.name = name
        self.idict = {}
        if odict != None:
            self.idict = dict((i,GeneSet(pair=(i,g))) for i,g in odict.items())
        elif gs != None:
            #FIXME check id if it already exists?
            self.idict = dict((g.id,g) for g in gs)

    def ids(self):
        return self.idict.keys()

    def get(self, a):
        return self.idict[a]

    def to_odict(self):
        """ Return gene sets in old dictionary format. """
        return dict(gs.to_odict() for gs in self.idict.values())

    def add(self, gs):
        if gs.id in idict:
            raise GeneSetIDException

    def __repr__(self):
        return "GeneSets(" + str(self.name) + ", " + str(self.idict) + ")"

def geneSetUnion():
    pass

def goGeneSets(org):
    """Returns gene sets from GO."""

    ontology = obiGO.Ontology.Load()
    annotations = obiGO.Annotations.Load(org, ontology=ontology)

    genesets = []

    for termn, term in ontology.terms.items():
        genes = annotations.GetAllGenes(termn)
        hier = ("GO", term.namespace)
        gs = GeneSet(id=termn, name=term.name, genes=genes, hierarchy=hier, organism=org) 
        genesets.append(gs)

    return GeneSets(gs=genesets)

def keggGeneSets(org):
    """
    Returns gene sets from KEGG pathways.
    """
    kegg = obiKEGG.KEGGOrganism(org)

    genesets = []
    for id in kegg.pathways():
        pway = obiKEGG.KEGGPathway(id)
        hier = ("KEGG",)
        gs = GeneSet(id=id, name=pway.title, genes=kegg.get_genes_by_pathway(id), hierarchy=hier, organism=org)
        genesets.append(gs)

    return GeneSets(gs=genesets)

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

def loadGMT(contents, name):
    """
    Eech line consists of tab separated elements. First is
    the geneset name, next is it's description. 
    
    For now the description is skipped.
    """
    def hline(s):
        tabs = [ strip(tab) for tab in s.split("\t") ]
        return  GeneSet(id=tabs[0], description=tabs[1], hierarchy=(name,), genes=tabs[2:])

    return GeneSets(gs=handleNELines(contents, hline))

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

def createCollection(lnf):
    """
    Input - list of tuples of geneset collection name and GMT
    file.
    """

    gen1 = {}

    for n,fn in lnf:

        if fn.lower()[-4:] == ".gmt": #format from webpage
            gen2 = loadGMT(open(fn,"rt").read(), fn).to_odict()
            gen1.update(gen2)

        elif n == fn and n[0] == ":":
            _,name,org = n.split(":")

            if name == "kegg":
                gen2 = keggGeneSets(org).to_odict()
                gen1.update(gen2)

            elif name == "go":
                gen2 = goGeneSets(org).to_odict()
                gen1.update(gen2)

            else:
                raise Exception("Wrong special name (%s)" % (name))            
            
        else:
            raise Exception("Can not recognize geneset (%s,%s)" % (n,fn))

    return gen1

def getCollectionFiles(path=collectionsPathname()):
    """ Get collections from a given path name. """

    def loadInfo(path):
        #TODO load info file
        return {}
        
    info = loadInfo(path)

    def okFile(fn):
        if fn.lower()[-4:] == ".gmt":
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
            if name.lower()[-4:] == ".gmt":
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
    #print keggGeneSets("sce").items()[:10]
    #col = collections([":go:sce"])
    #col = collections([":kegg:9606", ":go:9606"])
    #col = collections([":kegg:9606"])
    col = collections(["C2.CP.gmt"])
    #print col.items()[:10]
    import random
    print random.Random(0).sample(col.items(), 20)
