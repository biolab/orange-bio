
ncbi_geneinfo_tags = ("tax_id", "gene_id", "symbol", "locus_tag", "synonyms",
                      "dbXrefs", "chromosome", "map_location", "description", "type",
                      "symbol_from_nomenclature_authority", "full_name_from_nomenclature_authority",
                      "nomenclature_status", "other_designations", "modification_date")

ncbi_multiple_cardinality_tags = ("synonyms", "dbXrefs", "other_designations")

import os
import obiTaxonomy
import orngServerFiles

default_database_path = orngServerFiles.localpath("NCBI_geneinfo")

class GeneInfo(object):
    """ An object representing the NCBI information for a gene.
    """
    __slots__ = ncbi_geneinfo_tags
    def __init__(self, line):
        """ Construct the GeneInfo object from a line in the NCBI gene_info file
        """
        line = line.split("\t")
        for attr, value in zip(self.__slots__, line):
            if value == "-":
                value = None
            if attr in ncbi_multiple_cardinality_tags:
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

class NCBIGeneInfo(dict):
    _object_cache = {}
    def __init__(self, *args, **kwargs):
        """ An dictionary like object for accessing NCBI gene info
        Arguments::
                - *organsim*    Organism id

        Example::
            >>> info = NCBIGeneInfo("Homo sapiens")
        """
        if args and type(args[0]) in [str, unicode]:
            org = args[0]
            taxids = obiTaxonomy.to_taxid(org, mapTo=[obiTaxonomy.common_taxids()])
            if not taxids:
                taxids = set(obiTaxonomy.common_taxids()).intersection(obiTaxonomy.search(org))
            if len(taxids) == 0:
                raise obiTaxonomy.UnknownSpeciesIdentifier, org
            elif len(taxids) > 1:
                raise obiTaxonomy.MultipleSpeciesException, ", ".join(["%s: %s" % (id, obiTaxonomy.name(id)) for id in taxids])
            
            self.taxid = taxids.pop()
            if not os.path.exists(orngServerFiles.localpath("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)):
                orngServerFiles.download("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)
            file = open(orngServerFiles.localpath("NCBI_geneinfo", "gene_info.%s.db" % self.taxid), "rb")
            self.update(dict((line.split("\t", 3)[1], line) for line in file.read().split("\n") if line.strip() and not line.startswith("#")))

        else:
            dict.__init__(self, *args, **kwargs)
            
        # following is a temporary fix before gene name matcher is complete (then, this code is to be replaced)
        print self.keys()[:10]

        self.translate = dict([(self[k].symbol, k) for k in self.keys()])
        for k in self.keys():
            self.translate.update([(s, k) for s in self[k].synonyms if s not in self.translate] + \
                                  ([(self[k].locus_tag, k)] if self[k].locus_tag else [] ))

    @classmethod    
    def load(cls, file):
        """ A class method that loads gene info from file
        """
        if type(file) in [str, unicode]:
            file = open(file, "rb")
        return cls((line.split("\t", 3)[1], line) for line in file.read().split("\n") if line.strip() and not line.startswith("#"))
        
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
        id = self.translate.get(name, name)
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

    def search(self, string, exact=False):
        pass

    @staticmethod
    def get_geneinfo_from_ncbi(progressCallback=None):
        import urllib2, gzip
        from cStringIO import StringIO
        data = StringIO(urllib2.urlopen("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz").read())
        info = NCBIGeneInfo.load(gzip.GzipFile(None, "rb", fileobj=data))
        return info
        

 

"""
Gene matcher.

"Database" for each oranism is a list of sets of gene aliases.
"""

from collections import defaultdict
import os

BUFFER_PATH = None

def ignore_case(gs):
    """ Transform names in sets in list to lower case """
    return [ set([a.lower() for a in g]) for g in gs ]

def create_mapping(groups):
    """ 
    Returns mapping of aliases to the group index.

    Unpickling the results of this function (binary format)
    is slower than running it.
    """
    togroup = defaultdict(list)
    for i,group in enumerate(groups):
        for alias in group:
            togroup[alias].append(i)
    return togroup

def join_sets(set1, set2):
    """ 
    Joins two gene set mapping. 

    A group g1 from set1 is joined to a group of aliases g2 from set2, 
    if and only if there intersection between g1 and g2 is not empty. 
    Returned all joined groups + groups that were not matched (returned
    unchanged).

    The operation is commutatitve and associative.
    """

    cur = [ set(a) for a in set1 ]
    currentmap = create_mapping(cur)

    new = [] #new groups

    #remember used to find unused
    set1used = set() 
    set2used = set()

    for i, group in enumerate(set2):

        #find groups of aliases (from set1)  intersecting with a current
        #group from set2
        cross = set(reduce(lambda x,y: x+y, 
            [ currentmap[alias] for alias in group if alias in currentmap], []))

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
 
def join_sets_l(lsets):
    """
    Joins multiple gene set mappings. Since joining is associative, we
    can simply chain joining.
    """
    current = lsets[0]
    for b in lsets[1:]:
        current = join_sets(current, b)
    return current

class Matcher(object):

    ignore_case = True

    #def __init__(self, ignore_case=True):
    #    self.ignore_case = ignore_case

    def set_targets(self, tl):
        """
        Set list on gene names tl as targets of this gene matcher. 
        Abstract function.
        """
        notImplemented()

    def match(self, gene):
        """Returns matching target gene name."""
        notImplemented()

class MatcherSequence(Matcher):
    """
    Chaining of gene matchers.
    
    User defines the order of gene matchers. Each gene is goes through
    sequence of gene matchers until a match is found.
    """
    
    def __init__(self, matchers):
        self.matchers = matchers

    def match(self, gene):
        for matcher in self.matchers:
            m = matcher.matchOne(gene)
            if m != None:
                return m
        return None

    def set_targets(self, targets):
        for matcher in self.matchers:
            matcher.set_targets(targets)

def buffer_path():
    """ Returns buffer path. Ignore it optionally. """
    if BUFFER_PATH == None:
        import orngEnviron
        pth = os.path.join(orngEnviron.directoryNames["bufferDir"], 
            "gene_matcher")
        try:
            os.makedirs(pth)
        except:
            pass
        return pth
    else:
        return BUFFER_PATH


def auto_pickle(filename, version, func, *args, **kwargs):
    """
    Run function func with given arguments and save results to
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
            if versionF == version:
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
    Forges a new matcher from list of sets of given aliases.
    """
    def __init__(self, aliases):
        print "running parent constructior"
        self.aliases = aliases
        self.mdict = create_mapping(self.aliases)

    def to_ids(self, gene):
        if self.ignore_case:
            gene = gene.lower()
        return self.mdict[gene]

    def set_targets(self, targets):
        d = defaultdict(list)
        for target in targets:
            ids = self.to_ids(target)
            if ids != None:
                for id in ids:
                    d[id].append(target)
        self.to_targets = d

    def match(self, gene):
        inputgeneids = self.to_ids(gene)
        #return target genes with same ids
        return set( \
            reduce(lambda x,y:x+y, 
                [ self.to_targets[igid] for igid in inputgeneids ], [])) 


class MatcherAliasesPickled(MatcherAliases):
    """
    Gene matchers using pickling should extend this class.
    
    Loading is done in a lazy way. Therefore defining joined gene matchers
    does not force full loading of its component, if they are not needed
    (joined pickled file is already prepared).
    """
    
    def set_aliases(self, aliases):
        self.saved_aliases = aliases

    def get_aliases(self):
        if not self.saved_aliases: #loads aliases if not loaded
            self.aliases = self.load_aliases()
        return self.saved_aliases

    aliases = property(get_aliases, set_aliases)

    def get_mdict(self):
        """ Creates mdict. Aliases are loaded if needed. """
        if not self.saved_mdict:
            self.saved_mdict = create_mapping(self.aliases)
        return self.saved_mdict

    def set_mdict(self, mdict):
        self.saved_mdict = mdict

    mdict = property(get_mdict, set_mdict)

    def set_targets(self, targets):
        MatcherAliases.set_targets(self, targets)

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
        filename = os.path.join(buffer_path(), self.filename())
        return auto_pickle(filename, self.create_aliases_version(), 
            self.create_aliases)

    def __init__(self):
        self.aliases = []
        self.mdict = {}
        print self.filename() # test if valid filename can be built

class MatcherAliasesKEGG(MatcherAliasesPickled):

    def _organism_name(self, organism):
        """ Returns internal KEGG organism name. Used to define file name. """
        import obiKEGG #FIXME speed up name resolving
        org = obiKEGG.KEGGOrganism(organism)
        return org.org

    def create_aliases(self):
        organism = self._organism_name(self.organism)
        import obiKEGG
        org = obiKEGG.KEGGOrganism(self.organism)
        genes = org.api._genes[org.org]
        return ignore_case([ set([name]) | set(b.alt_names) for 
            name,b in genes.items() ])

    def create_aliases_version(self):
        return "1"

    def filename(self):
        return "kegg_" + self._organism_name(self.organism)

    def __init__(self, organism):
        self.organism = organism
        MatcherAliasesPickled.__init__(self)

class MatcherAliasesGO(MatcherAliasesPickled):

    def create_aliases(self):
        import obiGO
        annotations = obiGO.Annotations.Load(self.organism)
        names = annotations.geneNamesDict
        return ignore_case(map(set, list(set([ \
            tuple(sorted(set([name]) | set(genes))) \
            for name,genes in names.items() ]))))

    def filename(self):
        import obiGO #FIXME name resolving is too slow for now.
        ao = obiGO.Annotations.Load(self.organism)
        return "go_" + os.path.basename(ao.file)

    def create_aliases_version(self):
        return "1" #FIXME need support for GO versioning

    def __init__(self, organism):
        self.organism = organism
        MatcherAliasesPickled.__init__(self)

class MatcherAliasesPickledJoined(MatcherAliasesPickled):
    """
    Forges a new matcher by joining gene aliases from different sets.
    """

    def filename(self):
        filenames = [ mat.filename() for mat in self.matchers ]
        return "__".join(filenames)

    def create_aliases(self):
        return join_sets_l([ mat.aliases for mat in self.matchers ])

    def create_aliases_version(self):
        return [ mat.create_aliases_version() for mat in self.matchers ]

    def __init__(self, matchers):
        """ 
        Join matchers together. Groups of aliases are joined if
        they share a common name.
        """
        #FIXME: sorting of matchers to avoid multipying pickled files for
        #different orderings.
        self.matchers = matchers
        MatcherAliasesPickled.__init__(self)
        
if __name__ == '__main__':
    """
    info = NCBIGeneInfo("9606")
    gi = info(list(info)[0])
    print gi.tax_id, gi.synonyms, gi.dbXrefs, gi.symbol_from_nomenclature_authority, gi.full_name_from_nomenclature_authority
    """

    import obiGeneSets

    def testsets():
        return obiGeneSets.collections([":kegg:hsa", ":go:hsa"])

    def names1():
        import orange
        data = orange.ExampleTable("DLBCL.tab")
        return [ a.name for a in  data.domain.attributes ]
    
    genesets = auto_pickle("testcol", "3", testsets)
    names = auto_pickle("testnames", "4", names1)
    print names[:100]
 
    import time

    print "loading time needs to be decreased to minimum"
    t = time.time()
    mat = MatcherAliasesKEGG("human")
    print "kegg", time.time() - t
    t = time.time()
    mat2 = MatcherAliasesGO("goa_human")
    print "go", time.time() - t
    t = time.time()
    mat3 = MatcherAliasesPickledJoined([mat,mat2])
    print "join", time.time() - t

    print "using targets"

    mat.set_targets(names)
    mat2.set_targets(names)
    mat3.set_targets(names)

    import mMisc as m

    print "before genes"
    genes = set(m.flatten(genesets.values()[:100]))
    print len(genes)
    print "after genes"

    for g in sorted(genes):
        print "KEGG", g, mat.match(g)
        print "GO  ", g, mat2.match(g)
        print "JOIN", g, mat3.match(g)

