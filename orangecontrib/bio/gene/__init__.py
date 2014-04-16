from __future__ import absolute_import

import os, time

from Orange.orng import orngServerFiles

from .. import taxonomy as obiTaxonomy
from .. import kegg as obiKEGG

from .. import dicty as obiDicty
from .. import biomart as obiBioMart

from . import homology

default_database_path = orngServerFiles.localpath("NCBI_geneinfo")

class GeneInfo(object):
    """ An object representing the NCBI information for a gene.
    """

    NCBI_GENEINFO_TAGS = ("tax_id", "gene_id", "symbol", "locus_tag", "synonyms",
                          "dbXrefs", "chromosome", "map_location", "description", "type",
                          "symbol_from_nomenclature_authority", "full_name_from_nomenclature_authority",
                          "nomenclature_status", "other_designations", "modification_date")
    NCBI_MULTIPLE_CARDINALITY_TAGS = ("synonyms", "dbXrefs", "other_designations")
    
    __slots__ = NCBI_GENEINFO_TAGS
    def __init__(self, line):
        """ Construct the GeneInfo object from a line in the NCBI gene_info file
        """
        line = line.split("\t")
        for attr, value in zip(self.__slots__, line):
            if value == "-":
                value = None
            if attr in GeneInfo.NCBI_MULTIPLE_CARDINALITY_TAGS:
                value = value.split("|") if value != None else []
            setattr(self, attr, value)

    def __repr__(self):
        def format(value):
            if not value:
                return "-"
            elif type(value) == list:
                return "|".join(value)
            else:
                return value
        return "\t".join(format(getattr(self, slot)) for slot in self.__slots__)

    def __str__(self):
        return repr(self)

class GeneHistory(object):
    NCBI_GENE_HISTORY_TAGS = ("tax_id", "gene_id", "discontinued_gene_id", "discontinued_symbol", "discontinue_date")
    __slots__ = NCBI_GENE_HISTORY_TAGS
    def __init__(self, line):
        for attr, value in zip(self.__slots__, line.split("\t")):
            setattr(self, attr, value)
            
            
class NCBIGeneInfo(dict):
    TAX_MAP = {
            "2104": "272634",  # Mycoplasma pneumoniae
            "4530": "39947",  # Oryza sativa
            "5833": "36329",  # Plasmodium falciparum
            "4932": "559292",  # Saccharomyces cerevisiae
            }
       
    def __init__(self, organism, genematcher=None):
        """ An dictionary like object for accessing NCBI gene info
        Arguments::
                - *organism*    Organism id

        Example::
            >>> info = NCBIGeneInfo("Homo sapiens")
        """
        
        self.taxid = self.organism_name_search(organism)


        fname = orngServerFiles.localpath_download("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)
        file = open(fname, "rb")
        self.update(dict([(line.split("\t", 3)[1], line) for line in file.read().splitlines() if line.strip() and not line.startswith("#")]))

        # NOTE orig init time for gene matcher: 2.5s, new 4s: investigate the slowdown
        # NOTE matches are not the same because aliases are build a bit
        # differently (main name versus old aliases conflict!)

        self.matcher = genematcher
        if self.matcher == None:
            if self.taxid == '352472':
                self.matcher = matcher([GMNCBI(self.taxid), GMDicty(), [GMNCBI(self.taxid), GMDicty()]])
            else:
                self.matcher = matcher([GMNCBI(self.taxid)])

        #if this is done with a gene matcher, pool target names
        self.matcher.set_targets(self.keys())
        
    def history(self):
        if getattr(self, "_history", None) is None:
            fname = orngServerFiles.localpath_download("NCBI_geneinfo", "gene_history.%s.db" % self.taxid)
            try:
                self._history = dict([(line.split("\t")[2], GeneHistory(line)) for line in open(fname, "rb").read().splitlines()])
                
            except Exception, ex:
                print >> sys.srderr, "Loading NCBI gene history failed.", ex
                self._history = {}
        return self._history
        
    @classmethod
    def organism_version(cls, name):
        oname = cls.organism_name_search(name)
        #FIXME, dirty hack to ensure file id downloaded
        orngServerFiles.localpath_download("NCBI_geneinfo", "gene_info.%s.db" % oname) 
        return orngServerFiles.info("NCBI_geneinfo", "gene_info.%s.db" % oname)["datetime"]

    @classmethod
    def organism_name_search(cls, org):
        if org in cls.common_taxids():
            return org
        elif org in NCBIGeneInfo.TAX_MAP:
            return NCBIGeneInfo.TAX_MAP[org]

        taxids = obiTaxonomy.to_taxid(org, mapTo=cls.common_taxids())
        if not taxids:
            taxids = obiTaxonomy.search(org, onlySpecies=False)
            taxids = set(cls.common_taxids()).intersection(taxids) #onlySpecies=False needed to find correct dicty
        if len(taxids) == 0:
            raise obiTaxonomy.UnknownSpeciesIdentifier, org
        elif len(taxids) > 1:
            raise obiTaxonomy.MultipleSpeciesException, ", ".join(["%s: %s" % (id, obiTaxonomy.name(id)) for id in taxids])
        taxid = taxids.pop()
        return cls.TAX_MAP.get(taxid, taxid)

    @classmethod    
    def load(cls, file):
        """ A class method that loads gene info from file
        """
        if type(file) in [str, unicode]:
            file = open(file, "rb")
        return cls((line.split("\t", 3)[1], line) for line in file.read().splitlines() if line.strip() and not line.startswith("#"))
        
    def get_info(self, gene_id, def_=None):
        """ Search and return the GeneInfo object for gene_id
        """
        try:
            return self(gene_id)
        except KeyError:
            return def_
        
    def __call__(self, name):
        """ Search and return the GeneInfo object for gene_id
        """
        #id = self.translate.get(name, name)
        #print self.matcher.umatch(name), self.matcher.match(name)
        id = self.matcher.umatch(name)
        return self[id]

    def __getitem__(self, key):
#        return self.get(gene_id, self.matcher[gene_id])
        return GeneInfo(dict.__getitem__(self, key))

    def __setitem__(self, key, value):
        if type(value) == str:
            dict.__setitem__(self, key, value)
        else:
            dict.__setitem__(self, key, repr(value))

    def get(self, key, def_=None):
        try:
            return self[key]
        except KeyError:
            return def_

    def itervalues(self):
        for val in dict.itervalues(self):
            yield GeneInfo(val)

    def iteritems(self):
        for key, val in zip(self.iterkeys(), self.itervalues()):
            yield key, val

    def values(self):
        return list(self.itervalues())
    
    def items(self):
        return list(self.iteritems())

    @staticmethod
    def get_geneinfo_from_ncbi(file, progressCallback=None):
        import urllib2, gzip, shutil, tempfile
        from cStringIO import StringIO
        if isinstance(file, basestring):
            file = open(file, "wb")
        
        stream = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz")
        tmpfile = tempfile.TemporaryFile()
        shutil.copyfileobj(stream, tmpfile)
        tmpfile.seek(0)
        stream = gzip.GzipFile(None, "rb", fileobj=tmpfile)
        shutil.copyfileobj(stream, file)
        
    @staticmethod
    def get_gene_history_from_ncbi(file, progressCallback=None):
        import urllib2, gzip, shutil, tempfile
        from cStringIO import StringIO
        if isinstance(file, basestring):
            file = open(file, "wb")
        
        stream = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz")
        tmpfile = tempfile.TemporaryFile()
        shutil.copyfileobj(stream, tmpfile)
        tmpfile.seek(0)
        stream = gzip.GzipFile(None, "rb", fileobj=tmpfile)
        shutil.copyfileobj(stream, file)
        
    @classmethod
    def common_taxids(cls):
        taxids = obiTaxonomy.common_taxids()
        return [cls.TAX_MAP.get(id, id) for id in taxids if cls.TAX_MAP.get(id, id)]
    
    @classmethod
    def essential_taxids(cls):
        taxids = obiTaxonomy.essential_taxids()
        return [cls.TAX_MAP.get(id, id) for id in taxids if cls.TAX_MAP.get(id, id)]


class EnsembleGeneInfo(object):
#    BIO_MART_DATASETS = {"9606": "hsapiens"}
    DEF_ATTRS = ["ensembl_gene_id", "external_gene_id", "entrezgene"]
    BIOMART_CONF = {"9606": ("hsapiens_gene_ensembl", DEF_ATTRS) # + ["hgnc_symbol"])
                    }
    def __init__(self, organism, gene_matcher=None):
        self.organism = self.organism_name_search(organism) 
        self.load()
        if gene_matcher is None:
            self.gene_matcher = matcher([GMEnsembl(self.organism), GMNCBI(self.organism)])
        else:
            self.gene_matcher = gene_matcher
        
    @classmethod
    def organism_name_search(cls, name):
        taxids = obiTaxonomy.to_taxid(name, mapTo=cls.common_taxids())
        if len(taxids) == 1:
            return taxids[0]
        else:
            raise ValueError(name)
        
    @classmethod
    def common_taxids(cls):
        return ["3702", "9913", "6239", "7955", "9606", "7227",
                "10090", "10116", "4932", "4896", "8355"]
    
    @classmethod
    def version(cls):
        return "v1.0"
        
    def filename(self):
        return "ensembl_" + self.organism
    
    def create_info(self):
        if self.organism in self.BIOMART_CONF:
            dset_name, attrs = self.BIOMART_CONF[self.organism]
        else:
            dset_name, attrs = self.default_biomart_conf(self.organism)
        
        dataset = obiBioMart.BioMartDataset("ensembl", dset_name)
        table = dataset.get_example_table(attributes=attrs, unique=True)
        print len(table)
        table.save(dset_name + ".tab")
        from collections import defaultdict
        names = defaultdict(set)
        for ex in table:
            row = [str(v) for v in ex if not v.isSpecial() and str(v)]
            names[row[0]].update(row)
        return names
        
    def default_biomart_conf(self, taxid):
        name = obiTaxonomy.name(self.organism).lower()
        name1, name2 = name.split(" ")[: 2]
        dset_name = name1[0] + name2 + "_gene_ensembl"
        return dset_name, self.DEF_ATTRS
    
    def load(self):
        import cPickle
        dir = orngServerFiles.localpath("EnsembleGeneInfo")
        if not os.path.exists(dir):
            os.makedirs(dir)
        
        try:
            filename = orngServerFiles.localpath_download("EnsembleGeneInfo", self.filename())
            info = cPickle.load(open(filename, "rb"))
        except Exception, ex:    
            filename = orngServerFiles.localpath("EnsembleGeneInfo", self.filename())
            info = self.create_info()
            cPickle.dump(info, open(filename, "wb"))
            
        self.info = info
        
    def __getitem__(self, key):
        return self.info[key]
    def __contains__(self, key):
        return key in self.info
    def __len__(self):
        return len(self.info)
    def __iter__(self):
        return iter(self.info)
    def keys(self):
        return self.info.keys()
    def values(self):
        return self.info.values()
    def items(self):
        return self.info.items()
    def get(self, key, default=None):
        return self.info.get(key, default)
    
    def genes(self):
        return self.info.keys()
    
    def aliases(self):
        return [(key,) + tuple(value) for key, value in self.items()]
    
    def ensembl_id(self, name):
        return self.gene_matcher.umatch(name)
        
"""
Gene matcher.

"Database" for each organism is a list of sets of gene aliases.
"""

from collections import defaultdict
import os

gene_matcher_path = None

def ignore_case(gs):
    """ Transform names in sets in list to lower case """
    return [ set([a.lower() for a in g]) for g in gs ]

def create_mapping(groups, lower=False):
    """ 
    Returns mapping of aliases to the group index. If lower
    is True, lower case forms of gene aliases are mapped to indices.

    Unpickling the results of this function (binary format)
    is slower than running it.

    TIMING NOTES: 
    - lower() costs are neglible (< 10%)
    - building sets instead of lists also costs about 10 percent
    """
    togroup = defaultdict(set)

    # code duplicated because running a function in relatively expensive here.
    if lower: 
        for i,group in enumerate(groups):
            for alias in group:
                togroup[alias.lower()].add(i)
    else:
        for i,group in enumerate(groups):
            for alias in group:
                togroup[alias].add(i)

    return togroup

def join_sets(set1, set2, lower=False):
    """ 
    Joins two sets of gene set mappings. If lower is True, lower case
    forms of gene aliases are compared.

    A group g1 from set1 is joined to a group of aliases g2 from set2, 
    if the groups share at least one gene. 
    Returns all joined groups and groups that were not matched, which 
    remain unchanged.

    The operation both commutative and associative.
    """

    set1 = [ set(a) for a in set1 ]
    set2 = [ set(a) for a in set2 ]

    currentmap = create_mapping(set1, lower=lower)

    new = [] #new groups

    #remember used to find unused
    set1used = set() 
    set2used = set()

    fn = lambda x: x
    if lower:
        fn = lambda x: x.lower()

    for i, group in enumerate(set2):

        #find groups of aliases (from set1)  intersecting with a current
        #group from set2
        cross = reduce(set.union, 
            [ currentmap[fn(alias)] for alias in group if fn(alias) in currentmap], set())

        for c in cross:
            #print c, group & set1[c], group, set1[c]
            set1used.add(c)
            set2used.add(i)
            new.append(group | set1[c]) #add a union

    #add groups without matches (from both sets)
    set1add = set(range(len(set1))) - set1used
    set2add = set(range(len(set2))) - set2used

    for a in set1add:
        new.append(set1[a])
    for a in set2add:
        new.append(set2[a])

    return new
 
def join_sets_l(lsets, lower=False):
    """
    Joins multiple gene set mappings using join_sets function.
    """
    current = lsets[0]
    for b in lsets[1:]:
        current = join_sets(current, b, lower=lower)
    return current

class Matcher(object):
    """
    Matches an input gene to some target gene (set in advance).
    """

    def copy(self):
        return notImplemented()

    def __call__(self, targets):
        return self.set_targets(targets)

    def set_targets(self, targets):
        """ Set input list of gene names (a list of strings) as target genes.
        """
        notImplemented()

    def match(self, gene):
        """Return a list of target gene aliases which share a set of aliases with the input gene (can be empty)."""
        notImplemented()

    def umatch(self, gene):
        """Return a single (unique) matching target gene or None, if there are no matches or multiple matches."""
        mat = self.match(gene)
        return mat[0] if len(mat) == 1 else None

    def explain(self, gene):
        """ 
        Return gene matches with explanations as lists of tuples:
        a list of matched target genes and the corresponding set of gene aliases.
        """
        notImplemented()

def buffer_path():
    """ Returns buffer path from Orange's setting folder if not 
    defined differently (in gene_matcher_path). """
    if  gene_matcher_path == None:
        from Orange.orng import orngEnviron
        pth = os.path.join(orngEnviron.directoryNames["bufferDir"], 
            "gene_matcher")
        try:
            os.makedirs(pth)
        except:
            pass
        return pth
    else:
        return gene_matcher_path

def auto_pickle(filename, version, func, *args, **kwargs):
    """
    Run function func with given arguments and save the results to
    a file named filename. If results for a given filename AND
    version were already saved, just read and return them.
    """

    import cPickle as pickle

    output = None
    outputOk = False

    try:
        f = open(filename,'rb')

        try:
            versionF = pickle.load(f)
            if version == None or versionF == version:
                outputOk = True
                output = pickle.load(f)
        except:
            pass
        finally:
            f.close()

    except:
        pass

    if not outputOk:
        output = func(*args, **kwargs)

        #save output before returning
        f = open(filename,'wb')
        pickle.dump(version, f, -1)
        pickle.dump(output, f, -1)
        f.close()

    return output

class MatcherAliases(Matcher):
    """
    Genes matcher based on a list of sets of given aliases.

    Target genes belonging to same sets of aliases as the input gene are 
    returned as matching genes.

    """
    def __init__(self, aliases, ignore_case=True):
        self.aliases = aliases
        self.ignore_case = ignore_case
        self.mdict = create_mapping(self.aliases, self.ignore_case)

    def to_ids(self, gene):
        """ Return ids of sets of aliases the gene belongs to. """
        if self.ignore_case:
            gene = gene.lower()
        return self.mdict[gene]

    def set_targets(self, targets):
        """
        A reverse dictionary is made according to each target's membership
        in the sets of aliases.
        """
        d = defaultdict(list)
        #d = id: [ targets ], where id is index of the set of aliases
        for target in targets:
            ids = self.to_ids(target)
            if ids != None:
                for id in ids:
                    d[id].append(target)
        mo = MatchAliases(d, self)
        self.matcho = mo #backward compatibility - default match object
        return mo

    #this two functions are solely for backward compatibility
    def match(self, gene):
        return self.matcho.match(gene)
    def explain(self, gene):
        return self.matcho.explain(gene)

class Match(object):

    def umatch(self, gene):
        """Returns an unique (only one matching target) target or None"""
        mat = self.match(gene)
        return mat[0] if len(mat) == 1 else None
 
class MatchAliases(Match):

    def __init__(self, to_targets, parent):
        self.to_targets = to_targets
        self.parent = parent

    def match(self, gene):
        """
        The `gene` is first mapped to ids of sets of aliases which contain
        it. Target genes from the same sets of aliases are returned
        as input's match.
        """
        inputgeneids = self.parent.to_ids(gene)
        #return target genes with same ids
        return list(set( \
            reduce(lambda x,y: x+y, 
                [ self.to_targets[igid] for igid in inputgeneids ], [])))

    def explain(self, gene):
        inputgeneids = self.parent.to_ids(gene)
        return [ (self.to_targets[igid], self.parent.aliases[igid]) for igid in inputgeneids ]

class MatcherAliasesPickled(MatcherAliases):
    """
    Gene matchers based on sets of aliases supporting pickling should
    extend this class. Subclasses must define functions "filename", 
    "create_aliases_version" and "create_aliases". Those are crucial for
    pickling of gene aliases to work.

    Loading of gene aliases is done lazily: they are loaded when they are
    needed. Loading of aliases for components of joined matchers is often 
    unnecessary and is therefore avoided. 
    """
    
    def set_aliases(self, aliases):
        self.saved_aliases = aliases

    def get_aliases(self):
        if not self.saved_aliases: #loads aliases if not loaded
            self.aliases = self.load_aliases()
        #print "size of aliases ", len(self.saved_aliases)
        return self.saved_aliases

    aliases = property(get_aliases, set_aliases)

    def get_mdict(self):
        """ Creates mdict. Aliases are loaded if needed. """
        if not self.saved_mdict:
            self.saved_mdict = create_mapping(self.aliases, self.ignore_case)
        return self.saved_mdict

    def set_mdict(self, mdict):
        self.saved_mdict = mdict

    mdict = property(get_mdict, set_mdict)

    def set_targets(self, targets):
        return MatcherAliases.set_targets(self, targets)

    def filename(self):
        """ Returns file name for saving aliases. """
        notImplemented()
        
    def create_aliases_version(self):
        """ Returns version of the source database state. """
        notImplemented()

    def create_aliases(self):
        """ Returns gene aliases. """
        notImplemented()

    def load_aliases(self):
        fn = self.filename()
        ver = self.create_aliases_version() #if version == None ignore it
        if fn != None:
            if isinstance(fn, tuple): #if you pass tuple, look directly
               filename = fn[0]
            else:
               filename = os.path.join(buffer_path(), fn)
            return auto_pickle(filename, ver, self.create_aliases)
        else:
            #if either file version of version is None, do not pickle
            return self.create_aliases()

    def __init__(self, ignore_case=True):
        self.aliases = []
        self.mdict = {}
        self.ignore_case = ignore_case
        self.filename() # test if valid filename can be built


class MatcherAliasesKEGG(MatcherAliasesPickled):
    """ Alias: GMKEGG. 
    """

    def _organism_name(self, organism):
        return obiKEGG.organism_name_search(organism)

    def create_aliases(self):
        org = obiKEGG.KEGGOrganism(self.organism, genematcher=GMDirect())
        osets = org.gene_aliases()
        return osets

    def create_aliases_version(self):
        # use KEGG short release string (e.g. '66.0+')
        release = obiKEGG.KEGGOrganism.organism_version(self.organism) + ".2"
        release, _ = release.split("/")
        return release

    def filename(self):
        return "kegg_2_" + self._organism_name(self.organism)

    def __init__(self, organism, ignore_case=True):
        self.organism = organism
        MatcherAliasesPickled.__init__(self, ignore_case=ignore_case)


class MatcherAliasesFile(MatcherAliasesPickled):

    def create_aliases(self):
        canNotCreateButCanOnlyOpen()

    def create_aliases_version(self):
        return None

    def filename(self):
        return (self.filename_,)

    def __init__(self, filename, ignore_case=True):
        self.filename_ = filename
        MatcherAliasesPickled.__init__(self, ignore_case=ignore_case)


class MatcherAliasesGO(MatcherAliasesPickled):
    """ Alias: GMGO.
    """

    def _organism_name(self, organism):
        """ Returns internal GO organism name. Used to define file name. """
        from .. import obiGO
        return obiGO.organism_name_search(self.organism)

    def create_aliases(self):
        from .. import obiGO
        annotations = obiGO.Annotations(self.organism, genematcher=GMDirect())
        names = annotations.geneNamesDict
        return map(set, list(set([ \
            tuple(sorted(set([name]) | set(genes))) \
            for name,genes in names.items() ])))

    def filename(self):
        return "go_" + self._organism_name(self.organism)

    def create_aliases_version(self):
        from .. import obiGO
        return "v2." + obiGO.Annotations.organism_version(self.organism)

    def __init__(self, organism, ignore_case=True):
        self.organism = organism
        MatcherAliasesPickled.__init__(self, ignore_case=ignore_case)

class MatcherAliasesDictyBase(MatcherAliasesPickled):
    """ Alias: GMDicty.
    """

    def create_aliases(self):
        db = obiDicty.DictyBase()
        #db.info, db.mappings
        infoa = [ set([id,name]) | set(aliases) for id,(name,aliases,_) in db.info.items() ]
        mappingsa = [ set(filter(None, a)) for a in db.mappings ]
        joineda = join_sets(infoa, mappingsa, lower=True)
        return joineda

    def create_aliases_version(self):
        return "v1." + obiDicty.DictyBase.version()

    def filename(self):
        return "dictybase" 

    def __init__(self, ignore_case=True):
        MatcherAliasesPickled.__init__(self, ignore_case=ignore_case)

class MatcherAliasesNCBI(MatcherAliasesPickled):
    """ Alias: GMNCBI.
    """

    def _organism_name(self, organism):
        return NCBIGeneInfo.organism_name_search(organism)

    def create_aliases(self):
        ncbi = NCBIGeneInfo(self.organism, genematcher=GMDirect())
        out = []
        for k in ncbi.keys():
            out.append(set(filter(None, [k, ncbi[k].symbol, ncbi[k].locus_tag] + [ s for s in ncbi[k].synonyms ] )))
        return out

    def filename(self):
        return "ncbi_" + self._organism_name(self.organism)

    def create_aliases_version(self):
        return "v2." + NCBIGeneInfo.organism_version(self.organism)

    def __init__(self, organism, ignore_case=True):
        self.organism = organism
        MatcherAliasesPickled.__init__(self, ignore_case=ignore_case)
        
class MatcherAliasesAffy(MatcherAliasesPickled):
    def create_aliases(self):
        filename = orngServerFiles.localpath_download("Affy", self.organism + ".pickle")
        import cPickle
        return cPickle.load(open(filename, "rb"))
    
    def filename(self):
        return "affy_" + self.organism
    
    def create_aliases_version(self):
        orngServerFiles.localpath_download("Affy", self.organism + ".pickle")
        return orngServerFiles.info("Affy", self.organism + ".pickle")["datetime"]
        
    def __init__(self, organism, **kwargs):
        self.organism = organism
        MatcherAliasesPickled.__init__(self, **kwargs)
    
        
class MatcherAliasesEnsembl(MatcherAliasesPickled):
    """ A matcher for Ensemble ids. Alias: GMEnsemble.
    """
    DEF_ATTRS = ["ensembl_gene_id", "external_gene_id", "entrezgene"]
    # taxid: (dataset_name, [name_attr1, name_attr2 ...])
    BIOMART_CONF = {}
    def __init__(self, organism, **kwargs):
        self.organism = organism
        MatcherAliasesPickled.__init__(self, **kwargs)
        
    def filename(self):
        return "ensembl_" + self.organism
    
    def create_aliases_version(self):
        return "v1"
    
    def create_aliases(self):
        from .. import obiBioMart
        if self.organism in self.BIOMART_CONF:
            dset_name, attrs = self.BIOMART_CONF[self.organism]
        else:
            dset_name, attrs = self.default_biomart_conf(self.organism)
        
        dataset = obiBioMart.BioMartDataset("ensembl", dset_name)
        table = dataset.get_example_table(attributes=attrs)
        from collections import defaultdict
        names = defaultdict(set)
        for ex in table:
            row = [str(v) for v in ex if not v.isSpecial() and str(v)]
            names[row[0]].update(row)
        return names.values()
        
    def default_biomart_conf(self, taxid):
        name = obiTaxonomy.name(self.organism).lower()
        name1, name2 = name.split(" ")[: 2]
        dset_name = name1[0] + name2 + "_gene_ensembl"
        return dset_name, self.DEF_ATTRS
        

class MatcherAliasesPickledJoined(MatcherAliasesPickled):
    """
    Creates a new matcher by joining gene aliases from different data sets.
    Sets of aliases are joined if they contain common genes.

    The joined gene matcher can only be pickled if the source gene
    matchers are picklable.
    """

    def filename(self):
        # do not pickle if any is unpicklable
        try:
            filenames = [ mat.filename() for mat in self.matchers ]
            if self.ignore_case:
                filenames += [ "icj" ]
            return "__".join(filenames)
        except:
            return None

    def create_aliases(self):
        return join_sets_l([ mat.aliases for mat in self.matchers ], lower=self.ignore_case)

    def create_aliases_version(self):
        try:
            return "v4_" + "__".join([ mat.create_aliases_version() for mat in self.matchers ])
        except:
            return None

    def __init__(self, matchers):
        """ 
        Join matchers together. Groups of aliases are joined if
        they share a common name.

        If ignore_case is True, ignores case when joining gene aliases.
        """
        #FIXME: sorting of matchers to avoid multipying pickled files for
        #different orderings.
        self.matchers = matchers
        allic = set([ m.ignore_case for m in self.matchers ])
        if len(allic) > 1:
            notAllMatchersHaveEqualIgnoreCase()
        ignore_case = list(allic)[0]

        MatcherAliasesPickled.__init__(self, ignore_case=ignore_case)
        
class MatcherSequence(Matcher):
    """
    Each gene goes through sequence of gene matchers (in the same order
    as in the matchers arguments) until a match is found.
    """
    
    def __init__(self, matchers):
        self.matchers = matchers

    def set_targets(self, targets):
        ms = []
        targets = list(targets) #copy targets as multiple use would
                                #be problematic if a generator was passed
        for matcher in self.matchers:
            ms.append(matcher.set_targets(targets))
        om = MatchSequence(ms)
        self.matcho = om
        return om

    #this two functions are solely for backward compatibility
    def match(self, gene):
        return self.matcho.match(gene)

    def explain(self, gene):
        return self.matcho.explain(gene)

class MatchSequence(Match):

    def __init__(self, ms):
        self.ms = ms

    def match(self, gene):
        for match in self.ms:
            m = match.match(gene)
            if m: 
                return m
        return []

    def explain(self, gene):
        for match in self.ms:
            m = match.match(gene)
            if m: 
                return match.explain(gene)
        return []

class MatcherDirect(Matcher):
    """
    Directly match target names. Can ignore case. Alias: GMDirect.
    """

    def __init__(self, ignore_case=True):
        self.ignore_case = ignore_case

    def set_targets(self, targets):
        targets = list(targets) #would be problematic if generator was passed
                                #as it is used twice
        aliases = [ set([a]) for a in targets]
        self.am = MatcherAliases(aliases, ignore_case=self.ignore_case)
        self.matcho = self.am.set_targets(targets)
        return self.matcho

    #this two functions are solely for backward compatibility
    def match(self, gene):
        return self.matcho.match(gene)
    def explain(self, gene):
        return self.matcho.explain(gene)

               
GMDirect = MatcherDirect
GMKEGG = MatcherAliasesKEGG
GMGO = MatcherAliasesGO
GMNCBI = MatcherAliasesNCBI
GMDicty = MatcherAliasesDictyBase
GMAffy = MatcherAliasesAffy
GMEnsembl = MatcherAliasesEnsembl

def issequencens(x):
    return hasattr(x, '__getitem__') and not isinstance(x, basestring)

def matcher(matchers, direct=True, ignore_case=True):
    """
    Builds a new matcher from a list of gene matchers. Apply matchers in
    the input list successively until a match is found. If an element of
    `matchers` is a list, combine matchers in the sublist by joining overlapping 
    sets of aliases.

    :param list matchers: Gene matchers. 
    :param bool direct: If True, first try
      to match gene directly (a :obj:`MatcherDirect` is inserted in front of the
      gene matcher sequence).  
    :param bool ignore_case: passed to the added
      direct matcher.
    """
    seqmat = []
    if direct:
        seqmat.append(MatcherDirect(ignore_case=ignore_case))
    for mat in matchers:
        if issequencens(mat):
            mat = MatcherAliasesPickledJoined(list(mat))
            seqmat.append(mat)
        else:
            seqmat.append(mat)
    return MatcherSequence(seqmat)

if __name__ == '__main__':

    #m1 = matcher([[GMNCBI('44689'), GMDicty()]])
    #print m1.matchers[1].aliases[:100]

    #m2 = GMNCBI('Dictyostelium discoideum')
    #print m2.aliases



    """
    gi = info(list(info)[0])
    print gi.tax_id, gi.synonyms, gi.dbXrefs, gi.symbol_from_nomenclature_authority, gi.full_name_from_nomenclature_authority
    """

    #dobim z joinom prave stvari?

    import time
    from .. import obiGeneSets

    def testsets():
        return obiGeneSets.collections([":kegg:hsa", ":go:hsa"])

    def names1():
        import orange
        data = orange.ExampleTable("DLBCL.tab")
        return [ a.name for a in  data.domain.attributes ]

    def namesd():
        import orange
        data = orange.ExampleTable("dd_ge_biorep1.tab")
        names = [ ex["gene"].value for ex in data ]
        return names

    genesets = auto_pickle("testcol", "3", testsets)
    names = auto_pickle("testnames", "4", names1)
    names2 = auto_pickle("testnamesdicty", "4", namesd)

    info = NCBIGeneInfo('Dictyostelium discoideum')
    for a in names2:
        print a
        info.get_info(a)

    t = time.time()
    mat5 = matcher([[GMKEGG('human'),GMGO('human')]], direct=False, ignore_case=True)
    mat7 = GMDicty()
    mat8 = GMNCBI('Homo sapiens')
    print "initialized all", time.time()-t

    print "using targets"

    mat5.set_targets(names)
    mat7.set_targets(names)
    mat8.set_targets(names)

    print "before genes"
    genes = reduce(set.union, genesets.values()[:1000], set())
    genes = list(genes)
    print genes[:20]
    print len(genes)
    print "after genes"

    for g in sorted(genes):
        print "KGO ", g, mat5.match(g), mat5.explain(g)
        print "DICT", g, mat7.match(g)
        print "NCBI", g, mat8.match(g)


