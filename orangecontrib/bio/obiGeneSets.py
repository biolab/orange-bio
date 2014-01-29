"""
Maintainer: Marko Toplak
"""

from __future__ import absolute_import, with_statement

if __name__ == "__main__":
    __package__ = "Orange.bio"

import cPickle as pickle, os, tempfile, sys
from collections import defaultdict
import datetime

import Orange.core as orange
from Orange.orng import orngServerFiles

from . import obiGO, obiKEGG, obiTaxonomy

sfdomain = "gene_sets"

def nth(l,n):
    return [ a[n] for a in l]

from Orange.bio.geneset import GeneSet, GeneSets, GenesetRegException

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

def keggGeneSets(org):
    """
    Returns gene sets from KEGG pathways.
    """

    kegg = obiKEGG.KEGGOrganism(org)

    genesets = []
    for id in kegg.pathways():
        print id
        pway = obiKEGG.KEGGPathway(id)
        hier = ("KEGG","pathways")
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
    from . import obiDictyMutants
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
    from . import obiOMIM    
    genesets = [GeneSet(id=disease.id, name=disease.name, genes=obiOMIM.disease_genes(disease), hierarchy=("OMIM",), organism="9606",
                    link=("http://www.omim.org/entry/%s" % disease.id if disease.id else None)) \
                    for disease in obiOMIM.diseases()]
    return GeneSets(genesets)

def miRNAGeneSets(org):
    """
    Return gene sets from miRNA targets
    """
    from . import obimiRNA
    org_code = obiKEGG.from_taxid(org)
    link_fmt = "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s"
    mirnas = [(id, obimiRNA.get_info(id)) for id in obimiRNA.ids(org_code)]
    genesets = [GeneSet(id=mirna.matACC, name=mirna.matID, genes=mirna.targets.split(","), hierarchy=("miRNA", "Targets"),
                        organism=org, link=link_fmt % mirna.matID) for id, mirna in mirnas]
    return GeneSets(genesets)

def go_miRNASets(org, ontology=None, enrichment=True, pval=0.05, treshold=0.04):
    from . import obimiRNA, obiGO
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
        return GeneSet(id=tabs[0], description=tabs[1],
                                   hierarchy=("Custom Gene Sets",name), genes=tabs[2:])

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
    localfiles = set(orngServerFiles.listfiles(sfdomain))
    return [ filename_parse(fn) + \
        ((True,) if fn in localfiles else (False,)) for fn in gs_files ]

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

def list_all():
    """
    Return gene sets available in the local and ServerFiles repositories. 
    It returns a list of tuples of (hierarchy, organism, available_locally)
    """
    flist = list_local() + list_serverfiles()
    d = {}
    for h,o,local in flist:
        d[h,o] = min(local, d.get((h,o),True))
    return [ (h,o,local) for (h,o),local in d.items() ]

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
        taxname = obiTaxonomy.name(org)
        title = "Gene sets: " + ", ".join(hierarchy) + \
            ((" (" + taxname + ")") if org != None else "")
        tags = list(hierarchy) + [ "gene sets", taxname ] + obiTaxonomy.shortname(org) +\
            ([ "essential" ] if org in obiTaxonomy.essential_taxids() else [] )
        serverFiles.upload(sfdomain, fn, tfname, title, tags)
        serverFiles.unprotect(sfdomain, fn)
    finally:
        os.remove(tfname)

    update_server_list(serverFiles)

def register(genesets, serverFiles=None):
    """ Registers given :class:`GeneSets` locally.  The gene set is registered
    by the common hierarchy or organism (None if organisms are different).

    If :obj:`serverFiles` as a authenticated ServerFiles connection,
    the given gene sets are uploaded to the ServerFiles repository.  
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
            organismc = obiTaxonomy.to_taxid(organism)
            if len(organismc) == 1:
                organism = organismc.pop()
            else:
                exstr = "Could not interpret organism " + str(organism) + \
                      ". Possibilities: " + str(organismc) 
                raise NoGenesetsException(exstr)

    try:
        return load_local(hierarchy, organism)
    except NoGenesetsException:
        return load_serverfiles(hierarchy, organism)

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

if __name__ == "__main__":
    rsf = orngServerFiles.ServerFiles(username=sys.argv[1], password=sys.argv[2])
    upload_genesets(rsf)
    pass


