from __future__ import absolute_import, with_statement

import cPickle as pickle, os, tempfile, sys
from collections import defaultdict
import datetime

import Orange.core as orange
from Orange.orng import orngServerFiles

from .. import go as obiGO, kegg as obiKEGG, taxonomy as obiTaxonomy
from .. import dicty 
obiDictyMutants = dicty.phenotypes
from .. import omim as obiOMIM
from .. import go as obiGO

from . import transform

sfdomain = "gene_sets"

def nth(l,n):
    return [ a[n] for a in l]

class NoGenesetsException(Exception): pass

def goGeneSets(org):
    """Returns gene sets from GO."""
    ontology = obiGO.Ontology()
    annotations = obiGO.Annotations(org, ontology=ontology)

    genesets = []
    link_fmt = "http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s"
    for termn, term in ontology.terms.items():
        genes = annotations.GetAllGenes(termn)
        hier = ("GO", term.namespace)
        if len(genes) > 0:
            gs = GeneSet(id=termn, name=term.name, genes=genes, hierarchy=hier, organism=org, link=link_fmt % termn)
            genesets.append(gs)

    return GeneSets(genesets)

def strornone(x):
    return str(x) if x != None else x

def keggGeneSets(org):
    """
    Returns gene sets from KEGG pathways.
    """

    kegg = obiKEGG.KEGGOrganism(org)

    genesets = []
    for id in kegg.pathways():
        pway = obiKEGG.KEGGPathway(id)
        hier = ("KEGG","pathways")
        if pway.pathway_attributes():
            gs = GeneSet(id=id,
                                 name=pway.title,
                                 genes=kegg.get_genes_by_pathway(id),
                                 hierarchy=hier,
                                 organism=org,
                                 link=pway.link)
            genesets.append(gs)

    return GeneSets(genesets)

def dictyMutantSets():
    """
    Return dicty mutant phenotype gene sets from Dictybase
    """
    link_fmt = "http://dictybase.org/db/cgi-bin/dictyBase/SC/scsearch.pl?searchdb=strains&search_term=%s&column=all&B1=Submit" 
    #genesets = [GeneSet(id=mutant.name, name=mutant.descriptor, genes=obiDictyMutants.mutant_genes(mutant), hierarchy=("Dictybase", "Mutants"), organism="352472", # 352472 gathered from obiGO.py code_map -> Dicty identifier
    #                    link=(link_fmt % mutant.name if mutant.name else None)) \
    #                    for mutant in obiDictyMutants.mutants()]
 
    genesets = [GeneSet(id=phenotype, name=phenotype, genes=[obiDictyMutants.mutant_genes(mutant)[0] for mutant in mutants], hierarchy=("Dictybase", "Phenotypes"), organism="352472", # 352472 gathered from obiGO.py code_map -> Dicty identifier
                        link="") \
                        for phenotype, mutants in obiDictyMutants.phenotype_mutants().items()]

    return GeneSets(genesets)

def cytobandGeneSets():
    """
    Create cytoband gene sets from Stanford Microarray Database
    """
    import urllib2

    url = "http://www-stat.stanford.edu/~tibs/GSA/cytobands-stanford.gmt"
    stream = urllib2.urlopen(url)
    data = stream.read().splitlines()

    genesets = []
    for band in data:
        b = band.split("\t")
        genesets.append(GeneSet(id=b[0], name=b[1], genes=b[2:] if b[2:] else [], hierarchy=("Cytobands",), organism="9606", link=""))          

    return GeneSets(genesets)

def reactomePathwaysGeneSets():
    """
    Prepare human pathways gene sets from reactome pathways
    """
    import urllib
    import io
    from zipfile import ZipFile

    url = urllib.urlopen("http://www.reactome.org/download/current/ReactomePathways.gmt.zip")
    memfile = io.BytesIO(url.read())
    with ZipFile(memfile, "r") as myzip:
        f = myzip.open("ReactomePathways.gmt")
        content = f.read().splitlines()      

    genesets = [GeneSet(id=path.split("\t")[0], name=path.split("\t")[0], genes=path.split("\t")[2:] if path.split("\t")[2:] else [], hierarchy=("Reactome", "Pathways"), organism="9606", link="") for path in content]
    return GeneSets(genesets)


def omimGeneSets():
    """
    Return gene sets from OMIM (Online Mendelian Inheritance in Man) diseses
    """
    genesets = [GeneSet(id=disease.id, name=disease.name, genes=obiOMIM.disease_genes(disease), hierarchy=("OMIM",), organism="9606",
                    link=("http://www.omim.org/entry/%s" % disease.id if disease.id else None)) \
                    for disease in obiOMIM.diseases()]
    return GeneSets(genesets)

def miRNAGeneSets(org):
    """
    Return gene sets from miRNA targets
    """
    from .. import obimiRNA
    org_code = obiKEGG.from_taxid(org)
    link_fmt = "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s"
    mirnas = [(id, obimiRNA.get_info(id)) for id in obimiRNA.ids(org_code)]
    genesets = [GeneSet(id=mirna.matACC, name=mirna.matID, genes=mirna.targets.split(","), hierarchy=("miRNA", "Targets"),
                        organism=org, link=link_fmt % mirna.matID) for id, mirna in mirnas]
    return GeneSets(genesets)

def go_miRNASets(org, ontology=None, enrichment=True, pval=0.05, treshold=0.04):
    from .. import obimiRNA
    mirnas = obimiRNA.ids(int(org))
    if ontology is None:
        ontology = obiGO.Ontology()

    annotations = obiGO.Annotations(org, ontology=ontology)

    go_sets = obimiRNA.get_GO(mirnas, annotations, enrichment=enrichment, pval=pval, goSwitch=False)

    go_sets = obimiRNA.filter_GO(go_sets, annotations, treshold=treshold)

    link_fmt = "http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s"
    gsets = [GeneSet(id=key, name=ontology[key].name, genes=value, hierarchy=("miRNA", "go_sets",),
                        organism=org, link=link_fmt % key) for key, value in go_sets.items()]
    gset = GeneSets(gsets)
    return gset

import re

linkre = re.compile("(.*?)\s?(?:\[(https?://[^[^\]]*)\])?$")

def loadGMT(contents, name):
    """
    Eech line consists of tab separated elements. First is
    the geneset name, next is it's description.

    For now the description is skipped.

    Example Gene Set (.gmt) file:
    anti-liver_sw   anti-liver_sw   BMPR1A  APOD    WSB1    BMI1    SLC2A1  ...
    B-cells_sw  B-cells_sw  E2F5    NCF1    PALM2-AKAP2 IRF4    SLC2A1  ...
    Bladder_sw  Bladder_sw  PLCD4   ANGPTL1 LOC286191   ST0N1   LOC283904   ...
    cerebellum_sw   cerebellum_sw   C19orf43    LOC653464   KI110802    ...
    Cervix_sw   Cervix_sw   LAMA4   GSTM5   SNX19   DKK1    NT5E    ...
    """

    def hline(s):
        tabs = [tab.strip() for tab in s.split("\t")]
        groups = linkre.match(tabs[1]).groups()
        return GeneSet(id=tabs[0], description=groups[0], link=groups[1],
                                   hierarchy=("Custom",name), genes=tabs[2:])

    def handleNELines(s, fn):
        """
        Run function on nonempty lines of a string.
        Return a list of results for each line.
        """
        lines = (l.strip() for l in s.splitlines())
        return [fn(l) for l in lines if l]

    return GeneSets(handleNELines(contents, hline))

def getGenesetsStats(genesets):
    num_sets = len(genesets)
    unique_genes = len(set([gene for geneset in genesets for gene in geneset.genes]))
    genes_per_geneset = sum([len(geneset.genes) for geneset in genesets])/num_sets
    return num_sets, unique_genes, genes_per_geneset
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
    from Orange.orng import orngEnviron
    pth = os.path.join(orngEnviron.directoryNames["bufferDir"], "gene_sets_local")
    omakedirs(pth)
    return pth

def build_index(dir):
    """ Returns gene set availability index for some folder. """
    pass

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
    return [ filename_parse(fn) + (True,) for fn in gs_files ]

def remove_local(gene_set):
    """ Removes a given gene set from the local repository. """
    pth = local_path()
    gs_files = filter(is_genesets_file, os.listdir(pth)) 
    for setfile in gs_files:
        if setfile.__contains__(gene_set):
            setBgone = os.path.join(pth, setfile)
            os.remove(setBgone) 

def modification_date(file):
    t = os.path.getmtime(file)
    return datetime.datetime.fromtimestamp(t)

def list_serverfiles_from_flist(flist):
    gs_files = filter(is_genesets_file, flist)
    localfiles = os.listdir(orngServerFiles.localpath(sfdomain))
    localfiles = set(filter(is_genesets_file, localfiles))
    return [ filename_parse(fn) + \
        ((True,) if fn in localfiles else (False,)) for fn in set(gs_files) | localfiles ]

def list_serverfiles_conn(serverfiles=None):
    """ Returns available gene sets from the server files
    repository: a list of (hierarchy, organism, on_local) """
    if serverfiles == None:
        serverfiles = orngServerFiles.ServerFiles()
    flist = serverfiles.listfiles(sfdomain)
    return list_serverfiles_from_flist(flist)

def list_serverfiles():
    fname = orngServerFiles.localpath_download(sfdomain, "index.pck")
    flist = pickle.load(open(fname, 'r'))
    return list_serverfiles_from_flist(flist)

def list_all(org=None, local=None):
    """
    Return gene sets available in the local and ServerFiles repositories. 
    It returns a list of tuples of (hierarchy, organism, available_locally)

    Results can be filtered with the following parameters.

    :param str org: Organism tax id.
    :param bool local: Available locally.
    """
    flist = list_local() + list_serverfiles()
    d = {}
    for h,o,l in flist:
        d[h,o] = min(l, d.get((h,o),True))
    return [ (h,o,l) for (h,o),l in d.items() \
            if (local == None or l == local) and \
               (org == None or o == str(org))
        ]

def update_server_list(serverfiles_upload, serverfiles_list=None):
    if serverfiles_list == None:
        serverfiles_list = orngServerFiles.ServerFiles()

    flist = map(lambda x: filename(*x[:2]), list_serverfiles_conn(serverfiles_list))

    tfname = pickle_temp(flist)

    try:
        fn = "index.pck"
        title = "Gene sets: index"
        tags = [ "gene sets", "index", "essential" ]
        serverfiles_upload.upload(sfdomain, fn, tfname, title, tags)
        serverfiles_upload.unprotect(sfdomain, fn)
    except Exception,e:
        raise e
    finally:
        os.remove(tfname)

def _register_local(genesets):
    """ Registers using the common hierarchy and organism. """
    pth = local_path()
    org = genesets.common_org()
    hierarchy = genesets.common_hierarchy()
    fn = filename(hierarchy, org)

    with open(os.path.join(pth, fn), "w") as f:
        pickle.dump(genesets, f)

    return fn

def pickle_temp(obj):
    """ Pickle a file to a temporary file and returns its name """
    fd,tfname = tempfile.mkstemp()
    os.close(fd)
    f = open(tfname, 'wb')
    pickle.dump(obj, f)
    f.close()
    return tfname

def _register_serverfiles(genesets, serverFiles):
    """ Registers using the common hierarchy and organism. """
    org = genesets.common_org()
    hierarchy = genesets.common_hierarchy()
    fn = filename(hierarchy, org)

    #save to temporary file
    tfname = pickle_temp(genesets)

    try:
        if org != None:
            taxname = obiTaxonomy.name(org)
        title = "Gene sets: " + ", ".join(hierarchy) + \
            ((" (" + taxname + ")") if org != None else "")
        tags = list(hierarchy) + [ "gene sets" ] + ([ taxname ] if org != None else [])  + obiTaxonomy.shortname(org) +\
            ([ "essential" ] if org in obiTaxonomy.essential_taxids() else [] )
        serverFiles.upload(sfdomain, fn, tfname, title, tags)
        serverFiles.unprotect(sfdomain, fn)
    finally:
        os.remove(tfname)

    update_server_list(serverFiles)

def register(genesets, serverFiles=None):
    """ Registers given genesets locally.  The gene set is registered
    by the common hierarchy or organism (None if organisms are different).

    :param GeneSets genesets:
    :param serverFiles: If `serverFiles` is an authenticated ServerFiles connection,
        the input gene sets are uploaded to the repository.  
    """
    if serverFiles == None:
        _register_local(genesets)
    else:
        _register_serverfiles(genesets, serverFiles)

def build_hierarchy_dict(files):
    hierd = defaultdict(list)
    for ind,f in enumerate(files):
        hier, org = f
        for i in range(len(hier)+1):
            hierd[(hier[:i], org)].append(ind)
    return hierd

def load_local(hierarchy, organism):
    return load_fn(hierarchy, organism, list_local, 
        lambda h,o: os.path.join(local_path(), filename(h, o)))

def load_serverfiles(hierarchy, organism):
    return load_fn(hierarchy, organism, list_serverfiles, 
        lambda h,o: orngServerFiles.localpath_download(sfdomain, filename(h, o)))

def load_fn(hierarchy, organism, fnlist, fnget):
    files = map(lambda x: x[:2], fnlist())
    hierd = build_hierarchy_dict(files)
    out = GeneSets()
    matches = hierd[(hierarchy, organism)]
    if not matches:
        exstr = "No gene sets for " + str(hierarchy) + \
                " (org " + str(organism) + ")"
        raise NoGenesetsException(exstr)
    for (h, o) in [ files[i] for i in hierd[(hierarchy, organism)]]:
        fname = fnget(h, o)
        out.update(pickle.load(open(fname, 'r')))
    return out

def load(hierarchy, organism):
    """ First try to load from the local registered folder. If the file
    is not available, load it from the server files. """
    if organism != None:
        try:
            int(organism) #already a taxid
        except:
            organismc = obiTaxonomy.to_taxid(strornone(organism))
            if len(organismc) == 1:
                organism = organismc.pop()
            else:
                exstr = "Could not interpret organism " + str(organism) + \
                      ". Possibilities: " + str(organismc) 
                raise NoGenesetsException(exstr)

    try:
        return load_local(hierarchy, strornone(organism))
    except NoGenesetsException:
        return load_serverfiles(hierarchy, strornone(organism))

def collections(*args):
    """
    Load gene sets from various sources: GMT file, GO, KEGG, and others. 
    Return an instance of :class:`GeneSets`. 
    
    Each arguments specifies a gene set and can be either:
    
    * a filename of a GMT file,
    * a tuple (hierarchy, organism) (for example ``(("KEGG",), "10090")``), or
    * an instance of :class:`GeneSets`
    """
    result = GeneSets()

    for collection in args:
        try:
            result.update(collection)
        except (ValueError, TypeError):
            if issequencens(collection): #have a hierarchy, organism specification
                new = load(*collection)
                result.update(new)
            else:
                if collection.lower()[-4:] == ".gmt": #format from webpage
                    result.update(loadGMT(open(collection,"rt").read(), collection))
                else:
                    raise Exception("collection() accepts files in .gmt format only.")

    return result

def issequencens(x):
    "Is x a sequence and not string ? We say it is if it has a __getitem__ method and it is not an instance of basestring."
    return hasattr(x, '__getitem__') and not isinstance(x, basestring)

class TException(Exception): pass

def upload_genesets(rsf):
    """
    Builds the default gene sets and
    """
    orngServerFiles.update_local_files()

    genesetsfn = [ keggGeneSets, goGeneSets, miRNAGeneSets]
    organisms = obiTaxonomy.common_taxids()
    for fn in genesetsfn:
        for org in organisms:
            try:
                print "Uploading ORG", org, fn
                genesets = fn(org).split_by_hierarchy()
                for gs in genesets:
                    print "registering", gs.common_hierarchy()
                    register(gs, rsf) #server files
                    #register(gs)
                    print "successful", gs.common_hierarchy()
            except obiTaxonomy.UnknownSpeciesIdentifier:
                print "Organism ontology not available", org
            except GenesetRegException:
                print "Empty gene sets.", org


def only_option(a):
    if len(a) == 1:
        return list(a)[0]
    else:
        raise Exception()

class GenesetRegException(Exception): pass

class GeneSet(object):
    """ A single set of genes.
    """

    def __init__(self, genes=[], name=None, id=None, \
        description=None, link=None, organism=None, hierarchy=None, pair=None):
        """
        :param pair: Only for backward compatibility: convert a tuple (name, genes)
            into this object.
        """

        self.hierarchy = hierarchy     
        """ Hierarchy should be formated as a tuple, for example ``("GO", "biological_process")``"""

        self.genes = set(genes)
        """ A set of genes. Genes are strings. """

        self.name = name
        """ Gene set name. """

        self.id = id
        """ Short gene set ID. """

        self.description = description
        """ Gene set description. """

        self.link = link
        """ Link to further information about this gene set. """

        self.organism = organism
        """ Organism as a NCBI taxonomy ID. """

        if pair:
            self.id, self.genes = pair[0], set(pair[1])

    """
    the following functions are needed for sets of gene sets to be able
    to assess equality
    """

    def __hash__(self):
        return self.id.__hash__() + self.name.__hash__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def size(self):
        return len(self.genes)

    def cname(self, source=True, name=True):
        """ Return a gene set name with hieararchy. """
        oname = self.id
        if source and self.hierarchy:
            oname = "[ " + ", ".join(self.hierarchy) + " ] " + oname
        if name and self.name:
            oname = oname + " " + self.name
        return oname

    def to_odict(self, source=True, name=True):
        """
        For backward compatibility. Return a gene set as a tuple
        (id, list of genes).
        """
        return self.cname(source=source, name=name), self.genes

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

class GeneSets(set):
    """ A collection of gene sets: contains :class:`GeneSet` objects. 
    """
    
    def __init__(self, input=None):
        """
        If `input` is a dictionary, the gene sets are converted to the current format.
        """
        if input != None and len(input) > 0:
            self.update(input)

    def update(self, input):
        if input.__class__.__name__ == "GeneSets": #HACK: because location can appear different
            super(GeneSets, self).update(input)
        else:
            prepared_genesets = [] #parse them all before adding,
                                   #so that it fails on error
            if hasattr(input, "items"):
                for i, g in input.items():
                    prepared_genesets.append(GeneSet(pair=(i, g)))
            else:
                for i in input:
                    if isinstance(i, GeneSet):
                        prepared_genesets.append(i)
                    else:
                        i, g = i
                        prepared_genesets.append(GeneSet(pair=(i, g)))

            for g in prepared_genesets:
                self.add(g)

    def to_odict(self):
        """ Return gene sets in old dictionary format. """
        return dict(gs.to_odict() for gs in self)

    def set_hierarchy(self, hierarchy):
        """ Sets hierarchy for all gene sets. """
        for gs in self:
            gs.hierarchy = hierarchy

    def __repr__(self):
        return "GeneSets(" + set.__repr__(self) + ")"

    def common_org(self):
        """ Return a common organism. """
        if len(self) == 0:
            raise GenesetRegException("Empty gene sets.")

        organisms = set(a.organism for a in self)

        try:
            return only_option(organisms)
        except:
            raise GenesetRegException("multiple organisms: " + str(organisms))

    def hierarchies(self):
        """ Return all hierarchies. """
        if len(self) == 0:
            raise GenesetRegException("Empty gene sets.")
        return set(a.hierarchy for a in self)

    def common_hierarchy(self):
        """ Return a common hierarchy. """
        hierarchies = self.hierarchies()

        def common_hierarchy1(hierarchies):
            def hier(l): return set(map(lambda x: x[:currentl], hierarchies))
            currentl = max(map(len, hierarchies))
            while len(hier(currentl)) > 1:
                currentl -= 1
            return only_option(hier(currentl))

        return common_hierarchy1(hierarchies)

    def split_by_hierarchy(self):
        """ Split gene sets by hierarchies. Return a list of :class:`GeneSets` objects. """
        hd = dict((h,GeneSets()) for h in  self.hierarchies())
        for gs in self:
            hd[gs.hierarchy].add(gs)
        return hd.values()


if __name__ == "__main__":
    rsf = orngServerFiles.ServerFiles(username=sys.argv[1], password=sys.argv[2])
    upload_genesets(rsf)
    pass



