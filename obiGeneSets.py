"""
Getting genesets from KEGG and GO.

Maintainer: Marko Toplak
"""
from __future__ import with_statement

import obiKEGG, orange
import os
import obiGO
import cPickle as pickle
import orngServerFiles
import obiTaxonomy
import tempfile

sfdomain = "gene_sets"

def nth(l,n):
    return [ a[n] for a in l]

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

class GeneSets(dict):
    
    def __init__(self, odict=None, gs=None):
        """
        odict are genesets in old dict format.
        gs are genesets in new format
        """
        if odict != None:
            self.idict = dict((i,GeneSet(pair=(i,g))) for i,g in odict.items())
        elif gs != None:
            for g in gs:
                self[g.id] = g

    def to_odict(self):
        """ Return gene sets in old dictionary format. """
        return dict(gs.to_odict() for gs in self.values())

    def set_hierarchy(self, hierarchy):
        """ Sets hierarchy for all gene sets """
        for gs in self.values():
            gs.hierarchy = hierarchy

    def __repr__(self):
        return "GeneSets(" + str(self.name) + ", " + dict.__repr__(self) + ")"

    def common_org(self):
        if len(self) == 0:
            raise GenesetRegException("empty gene set")

        organisms = set(a.organism for a in self.values())

        try:
            return only_option(organisms)
        except:
            raise GenesetRegException("multiple organisms: " + str(organisms))

    def common_hierarchy(self):
        if len(self) == 0:
            raise GenesetRegException("empty gene set")

        hierarchies = set(a.hierarchy for a in self.values())

        def common_hierarchy1(hierarchies):
            def hier(l): return set(map(lambda x: x[:currentl], hierarchies))
            currentl = max(map(len, hierarchies))
            while len(hier(currentl)) > 1:
                currentl -= 1
            return only_option(hier(currentl))

        return common_hierarchy1(hierarchies)

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

def loadGMT(contents, name):
    """
    Eech line consists of tab separated elements. First is
    the geneset name, next is it's description. 
    
    For now the description is skipped.
    """
    def hline(s):
        tabs = [ tab.strip() for tab in s.split("\t") ]
        return  GeneSet(id=tabs[0], description=tabs[1], hierarchy=(name,), genes=tabs[2:])

    def handleNELines(s, fn):
        """
        Run function on nonempty lines of a string.
        Return a list of results for each line.
        """
        lines = s.split("\n")
        lines = [ l.strip() for l in lines ]
        lines = filter(lambda x: x != "", lines)
        return [ fn(l) for l in lines ]

    return GeneSets(gs=handleNELines(contents, hline))

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

def getCollectionFiles(path=None):
    """ Get collections from a given path name. """

    if path == None:
        path = local_path()

    def okFile(fn):
        if fn.lower()[-4:] == ".gmt":
            return True
        return False

    files = sorted(filter(okFile, os.listdir(path)))

    #print files, os.listdir(path)

    out = []

    for file in files:
        fn = os.path.join(path, file)
        name = file
        if name == file:
            #remove suffix if no info
            if name.lower()[-4:] == ".gmt":
                name = name[:-4]

        out.append( (name, fn) )

    return out

"""
We have multiple paths for gene set data:
buffer/bigfiles/gene_sets
and
buffer/gene_sets_local
both have available.txt
"""

def omakedirs(dir):
    try:
        os.makedirs(dir)
    except OSError:
        pass

def local_path():
    """ Returns local path for gene sets. Creates it if it does not exists
    yet. """
    import orngEnviron
    pth = os.path.join(orngEnviron.directoryNames["bufferDir"], "gene_set_local")
    omakedirs(pth)
    return pth

def build_index(dir):
    """ Returns gene set availability index for some folder. """
    pass

class GenesetRegException(Exception): pass

def only_option(a):
    if len(a) == 1:
        return list(a)[0]
    else:
        raise Exception()

def filename(hierarchy, organism):
    """ Obtain a filename for given hierarchy and organism. """
    return "gs_" + "_._".join(hierarchy + \
        (organism if organism != None else "",)) + ".pck"

def filename_parse(fn):
    """ Returns a hierarchy and the organism from the filename."""
    fn = fn[3:-4]
    parts = fn.split("_._")
    hierarchy = tuple(parts[:-1])
    org = parts[-1] if parts[-1] != "" else None
    return hierarchy, org

def is_genesets_file(fn):
    return fn.startswith("gs_") and fn.endswith(".pck")

def list_local():
    """ Returns available gene sets from the local repository:
    a list of (hierarchy, organism, on_local) """
    pth = local_path()
    gs_files = filter(is_genesets_file, os.listdir(pth))
    return [ filename_parse(fn) + (1,) for fn in gs_files ]
    
def list_serverfiles_conn(serverfiles=None):
    """ Returns available gene sets from the server files
    repository: a list of (hierarchy, organism, on_local) """
    if serverfiles == None:
        serverfiles = orngServerFiles.ServerFiles()
    return filter(is_genesets_file,
        serverfiles.listfiles(sfdomain))

def update_server_list(serverfiles_upload, serverfiles_list=None):
    if serverfiles_list == None:
        serverfiles_list = orngServerFiles.ServerFiles()
    flist = list_serverfiles_conn(serverfiles_list)

    tfname = pickle_temp(flist)
    
    try:
        fn = "index.pck"
        title = "Gene sets: index"
        tags = [ "gene sets", "index", "essential" ]
        serverfiles_upload.upload(sfdomain, fn, tfname, title, tags)
        serverfiles_upload.unprotect(sfdomain, fn)
    except Exception as e:
        raise e
    finally:
        os.remove(tfname)

def list_serverfiles():
    fname = orngServerFiles.localpath_download(sfdomain, "index.pck")
    return pickle.load(open(fname, 'r'))

def register_local(genesets):
    """ Registers using the common hierarchy and organism. """
    pth = local_path()

    org = genesets.common_org()
    hierarchy = genesets.common_hierarchy()

    fn = filename(hierarchy, org)
    with open(os.path.join(pth, fn), "w") as f:
        pickle.dump(genesets, f)

    return fn

def pickle_temp(obj):
    """ Pickle a file to a temporary file returns its name """
    fd,tfname = tempfile.mkstemp()
    os.close(fd)
    f = open(tfname, 'wb')
    pickle.dump(obj, f)
    f.close()
    return tfname

def register_serverfiles(genesets, serverFiles):
    """ Registers using the common hierarchy and organism. """
    org = genesets.common_org()
    hierarchy = genesets.common_hierarchy()
    fn = filename(hierarchy, org)

    #save to temporary file
    tfname = pickle_temp(genesets)
    
    try:
        title = "Gene sets: " + ", ".join(hierarchy) + \
            ((" (" + obiTaxonomy.name(org) + ")") if org != None else "")
        tags = list(hierarchy) + [ "gene sets", obiTaxonomy.name(org) ]
        serverFiles.upload(sfdomain, fn, tfname, title, tags)
        serverFiles.unprotect(sfdomain, fn)
    except Exception as e:
        raise e
    finally:
        os.remove(tfname)

    update_server_list(serverFiles)

def register(genesets, serverFiles=None):
    """
    Hierarchy is induced from the gene set names.
    """
    if serverFiles == None:
        register_local(genesets)
    else:
        register_serverfiles(genesets, serverFiles)

def collections(l=[], default=False, path=None):
    """
    Input is a list of collections.
    Default - if default collections are included.
    Input is a list of names. If names match to any names in path,
    they are taken. If not, file with that name is regarded as
    a filename of gene set colections
    """
    collections = getCollectionFiles(path)

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
    #col = collections(["C2.CP.gmt"])
    #print col.items()[:10]
    import sys
    gs = goGeneSets("9606")
    register_local(gs)
    rsf = orngServerFiles.ServerFiles(username=sys.argv[1], password=sys.argv[2])
    register_serverfiles(gs, rsf)
    print list_serverfiles_conn()
    print "Server list from index", list_serverfiles()
