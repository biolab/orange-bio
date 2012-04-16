"""
=======
obiKEGG
=======
obiKEGG module is an interface to `Kyoto Encyclopedia of Genes and Genomes
<http://www.genome.jp/kegg/>`_ that allows easy access to KEGG pathway
and genes data.

Example ::
    
    >>> genes = KEGGOrganism("human")
    >>> genes
"""

from __future__ import absolute_import

try:
    import Image, ImageDraw, ImageMath
except:
    pass

import cStringIO
import math
import time
import os, sys, tarfile
import re

from cPickle import load, loads, dump
from collections import defaultdict

import xml.dom.minidom as minidom

from functools import partial, wraps
import urllib2
import cPickle

from Orange.orng import orngEnviron, orngServerFiles

from . import obiData, obiProb, obiTaxonomy

DEFAULT_DATABASE_PATH = orngServerFiles.localpath("KEGG")
KEGG_FTP_PATH = "ftp://ftp.genome.jp/pub/kegg/"

_join = lambda *args: os.path.join(DEFAULT_DATABASE_PATH, *args)

def special_tags(domain, filename):
    return dict([tuple(tag.split(":", 1)) for tag in orngServerFiles.info(domain, filename)["tags"] \
                 if tag.startswith("#") and ":" in tag])
def loads(domain, filename, version=None):
    """ A function decorator for a function that loads data from server files. 
    Makes sure the file is downloaded and of the right version
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                ## First try downloading the required server file  
                file = filename if isinstance(filename, basestring) else filename(*args, **kwargs)
                orngServerFiles.localpath_download(domain, file)
                if version is not None and special_tags(domain, file).get("#version", str(version)) < str(version):
                    orngServerFiles.update(domain, file) 
            except Exception, ex:
                print ex, args, kwargs
            return func(*args, **kwargs)
        return wrapper
    return decorator
            
def deprecated(first_arg, text=None):
    import warnings
    def wrapper(*args, **kwargs):
        warnings.warn_explicit(("Call to deprecated function %s. " % first_arg.__name__) + (text or ""), 
                               category=DeprecationWarning,
                               filename=first_arg.func_code.co_filename,
                               lineno=first_arg.func_code.co_firstlineno + 1
                      )
#        warnings.warn(("Call to deprecated function %s." % first_arg.__name__) + (text or ""), 
#                               category=DeprecationWarning,;
#                               stacklevel =  -1 #-2 if text is None else -3
#                        )
        return first_arg(*args, **kwargs)
    if isinstance(first_arg, basestring): #deprecated called with a string as first argument
        return partial(deprecated, text=first_arg)
    else:
        return wraps(first_arg)(wrapper)
        
def cached(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        sig = args + tuple(sorted(kwargs.items()))
        if sig in wrapper:
            pass
    return wrapper
            
def cached_property(func):
    @property
    @wraps(func)
    def wrapper_func(self):
        cache_name = "_" + func.__name__
        if not hasattr(self, cache_name):
            setattr(self, cache_name, func(self))
        return getattr(self, cache_name)
    return wrapper_func

def cached_method(func, cache_name="_cached_method_cache", store=None):
    def wrapper(self, *args, **kwargs):
        sig = (func.__name__,) + args + tuple(sorted(kwargs.items()))
        if not hasattr(self, cache_name):
            setattr(self, cache_name, store() if store is not None else {})
        if sig not in getattr(self, cache_name):
            getattr(self, cache_name)[sig] = func(self, *args, **kwargs)
        return getattr(self, cache_name)[sig]
    return wrapper

class pickle_cache(dict):
    def __init__(self, filename, sync=True):
        dict.__init__(self)
        self.filename = filename
        if os.path.exists(self.filename):
            try:
                self.update(cPickle.load(open(self.filename, "rb")))
            except Exception, ex:
                print self.filename, "loading failed", ex
                
    def __setitem__(self, key, value):
        super(pickle_cache, self).__setitem__(key, value)
        self.sync()
        
    def __delitem__(self, key):
        super(pickle_cache, self).__delitem__(key)
        self.sync()
        
    def sync(self):
#        print "sync", self, self.filename
        cPickle.dump(dict(self.items()), open(self.filename, "wb"), 2)
        
def persistent_cached_method(func):
    return cached_method(func, cache_name="_persistent_cached_method_cache")

def persistent_cached_class(cls, store=pickle_cache):
    def cache_loader(self):
        if not hasattr(self, "_cache_filename"):
            setattr(self, "_cache_filename", getattr(self, "filename") + ".cache")
#        print "creating cache", self._cache_filename 
        try:
#            print os.stat(self.filename).st_mtime, os.stat(self._cache_filename).st_mtime
            if os.stat(self.filename).st_mtime > os.stat(self._cache_filename).st_mtime:
#                print "cache obsolete"
                os.remove(self._cache_filename)
        except Exception, ex:
            pass
        return store(self._cache_filename)
    cls._persistent_cached_method_cache = cached_property(cache_loader)
    return cls

class cache_manager(object):
    """ takes an fuction that returns a dict like object (cache)
    """
    def __init__(self, opener):
        self.opener = opener
        self.name = opener.__name__
        
    def __call__(self, *args, **kwargs):
        return self.opener(self, *args, **kwargs)
            
    def __get__(self, instance, type_=None):
        """ Return an cache store from instance """
        ## TODO: check if instance is a type
        if getattr(instance, self.name):
            return getattr(instance, self.name)
        return self.opener.__get__(instance, type)
    
    def cached_method(manager, wrapped):
        def wrapper(self, *args, **kwargs):
            sig = (wrapped.__name__,) + args + tuple(sorted(kwargs.items()))
            cache = getattr(self, manager.__name__)
            if sig not in cache:
                cache[sig] = func(self, *args, **kwargs)
            return cache[sig]

def proxy_dict_decorator(cls, proxy_name):
    for name in ["__len__", "__contains__", "__getitem__", "__iter__", "get", "has_key", "items", 
                 "iteritems", "iterkeys", "itervalues", "keys", "values", "update"]:
        def proxy(func):
            def f(self, *args, **kwargs):
                return func(getattr(self, proxy_name), *args, **kwargs)
            return f
        setattr(cls, name, proxy(getattr(dict, name)))
    return cls

def _create_path_for_file(target): 
    try:
        os.makedirs(os.path.dirname(target))
    except OSError:
        pass

def downloader(func):
    def wrapper(*args, **kwargs):
        import shutil
        def download(src, dst):
            if isinstance(dst, basestring):
                _create_path_for_file(dst)
                dst = open(dst, "wb")
            shutil.copyfileobj(src, dst)
        jobs = func(*args, **kwargs)
        for src, dst in [jobs] if type(jobs) == tuple else jobs:
            download(src, dst)
    return wrapper

from Orange.orng.orngMisc import with_gc_disabled

SLINE, MLINE, LINELIST = range(1, 4)

class _field_mutators(object):
    @classmethod
    def deindent(cls, text, field):
        pad = len(text) - len(text[len(field):].lstrip())
        return "\n".join([line[pad:] for line in text.splitlines()])
    
    @classmethod
    def to_list(cls, text, field):
        return [line for line in text.splitlines() if line.strip()]
    
    @classmethod
    def to_id_list(cls, text, field):
        return [tuple(line.split(": ", 1)) for line in cls.to_list(cls.deindent(text, field), field)]
    
    @classmethod
    def to_dict(cls, text, field):
        text = cls.deindent(text, field).replace("\n ", " ")
        lines = [t.split(": ",1) for t in text.split("\n")]
        return dict([(key, [t.strip() for t in value.split(" ") if t]) for key, value in lines])
    
    @classmethod
    def to_ids(cls, text, field):
        lines = [line.split(" ", 2)[:2] for line in cls.deindent(text, field).split("\n") if ":" in line.split(" ", 1)[0]]
        return [db.lower() + id for db, id in lines] 
#        ids = reduce(lambda col, line:1, text.split("\n"), [])
        return ids
    default = deindent
#    entry = classmethod(lambda cls, text, field: cls.deindent(text, field).split()[0])
    orthology = to_id_list
    pathway = classmethod(lambda cls, text, field: [t[1].split()[0] for t in cls.to_id_list(text, field)])
    codon_usage = classmethod(lambda cls, text, field: text.replace(field, " " * len(field)))
    structure = to_id_list
    motif = to_dict
    pathway = to_ids
    orthology = to_ids
    
    @classmethod
    def genes(cls, text, field):
        genes = [line.split(": ", 1) for line in cls.deindent(text, field).replace("\n    ", "").splitlines()]
        return dict([(org.lower(), [gene.split("(")[0] for gene in genes.split(" ")]) for org, genes in genes])
    
    @classmethod
    def dblinks(cls, text, field):
        links = [line.split(": ", 1) for line in cls.deindent(text, field).replace("\n   ", "").splitlines()]
        return dict([(db.lower(), [name for name in links.split(" ")]) for db, links in links])
    
    dblinks = to_dict

def entry_decorate(cls):
    reserved_map = {"class": "class_"}
    for field in cls.FIELDS:
        mutator = getattr(_field_mutators, field.lower(), _field_mutators.default)
        def construct(field, mutator):
            def get(self):
                text = self.get(field)
                return mutator(text, field) if text else None
            return get
        if not hasattr(cls, field.lower()):
            setattr(cls, reserved_map.get(field.lower(), field.lower()), property(construct(field, mutator)))
        else:
            import warnings
            warnings.warn(str(cls)+ "already contains " + field.lower())
    return cls

def borg_class(cls):
    borg__init__ = cls.__init__
    cls._borg_instance = None
    def __init__(self, *args, **kwargs):
        if not cls._borg_instance:
            borg__init__(self, *args, **kwargs)
            cls._borg_instance = args + tuple(sorted(kwargs.items())), self
        elif cls._borg_instance[0] == args + tuple(sorted(kwargs.items())):
            self.__dict__ = cls._borg_instance[1].__dict__
        else:
            print "Warning. Atempting to create an individual instance of '%s'.\nWe are borg. You will be assimilated. Resistance is futile." % cls.__name__
            self.__dict__ = cls._borg_instance[1].__dict__
    cls.__init__ = __init__
    return cls

def defaultpath(func):
#    @wraps
    def wrapper(*args, **kwargs):
        if "path" not in kwargs:
            kwargs["path"] = DEFAULT_DATABASE_PATH
        return func(*args, **kwargs)
    return wrapper
             
        
class KEGGDBEntry(object):
    FIELDS = []
    MULTIPLE_FIELDS = []
    def __init__(self, entrytext, index=None):
        self.entrytext = entrytext
        if not self.FIELDS:
            self.FIELDS = set([line.split()[0] for line in entrytext.split("\n") if line.strip() and not line.startswith(" ")])
        self._index(index)
        
    def _indices(self, field):
        text = self.entrytext
        fieldlen = len(field)
        return reduce(lambda indices,s: indices.append(text.find(field, indices[-1] + fieldlen)) or indices, range(text.count(field)), [-fieldlen])[1:]
         
    def _index(self, index=None):
        if index is not None:
            index, self.partitions = index
            self.index = dict(zip(self.FIELDS, index))
            return 
        index = dict([(field, self._indices(field)) for field in self.FIELDS])
        sorted_indices = sorted(reduce(set.union, index.values(), set())) + [-1]
        self.partitions = sorted_indices #[(sorted_indices[i], sorted_indices[i+1]) for i in range(len(sorted_indices) - 1)]
        self.index = dict([(key, [sorted_indices.index(ind) for ind in indices]) for key, indices in index.items()])
        
    def section_iter(self, field):
        for ind in self.index[field]:
            yield self.entrytext[self.partitions[ind]: self.partitions[ind+1]]
        
    def get(self, key, default=None):
        if key not in self.FIELDS:
            return default
        sections = [section for section in self.section_iter(key)]
        if key in self.MULTIPLE_FIELDS:
            return sections
        elif sections:
            return sections[0]
        else: 
            return default
        
    def entry_key(self):
        return self.entry.split()[0]
        
class KEGGDataBase(object):
    Entry = KEGGDBEntry
    VERSION = "v1.0"
    def __init__(self, file=None):
        self.load(file)
        
    @with_gc_disabled
    def load(self, file=None):
        if file is None:
            file = self.FILENAME % {"path":DEFAULT_DATABASE_PATH}
        self.filename = file
        data = open(file, "rb").read().split("///\n")
        file_timestamp = os.stat(self.filename).st_mtime
#        print file_timestamp
        build_index = False
#        if os.path.exists(self.filename + ".index"):
        try:
            version, timestamp, index = cPickle.load(open(self.filename + ".index", "rb"))
            assert(version == self.VERSION and timestamp == file_timestamp)
        except Exception, ex:
#            print "error", ex
            build_index = True
            index = [None] * len(data)
        self.entrys = [self.Entry(entry, ind) for entry, ind in zip(data, index) if entry.strip()]
        if build_index:
#            self.entrys = [self.Entry(entry) for entry in data if entry.strip()]
            index = [([e.index[field] for field in e.FIELDS], e.partitions) for e in self.entrys]
            cPickle.dump((self.VERSION, file_timestamp, index), open(self.filename + ".index", "wb"), cPickle.HIGHEST_PROTOCOL)
        
        self.entry_dict = dict([(entry.entry_key(), entry) for entry in self.entrys])
        
proxy_dict_decorator(KEGGDataBase, "entry_dict")
        
class KEGGGeneEntry(KEGGDBEntry):
    FIELDS = ["ENTRY", "NAME", "DEFINITION","ORTHOLOGY", "PATHWAY", "CLASS", "POSITION",
              "MOTIF", "DBLINKS", "STRUCTURE", "CODON_USAGE", "AASEQ", "NTSEQ"]
    
    def aliases(self):
        return [self.entry_key()] + (self.name.split(",") if self.name else []) + [link[1][0] for link in self.dblinks.items() if self.dblinks]
    
    @property
    def alt_names(self):
        return self.aliases()
        
entry_decorate(KEGGGeneEntry)
    
class KEGGGenes(KEGGDataBase):
    Entry = KEGGGeneEntry
    FILENAME = "genes/organisms/%(org_code)s/%(org_name)s.ent"
    
    @loads("KEGG", lambda self, org: "kegg_genes_%s.tar.gz" % org)
    def load(self, org):
        super(KEGGGenes, self).load(_join(self.filename(org)))
        self.entry_dict = dict([(org + ":" + key, value) for key, value in self.entry_dict.items()])
        self.org_code = org
        
    def gene_aliases(self):
        aliases = {}
        for entry in self.entrys:
            aliases.update(dict.fromkeys(entry.aliases(), self.org_code + ":" + entry.entry_key()))
        return aliases
            
    @classmethod
    def filename(cls, org):
        return cls.FILENAME % {"org_code":org, "org_name":KEGGGenome()[org].name.split(",")[1].strip()}
     
    @classmethod
    @downloader
    def download(cls, org):
#        local_dir = orngServerFiles.localpath("KEGG") if local_dir is None else local_dir
        filename = cls.filename(org)
        return (urllib2.urlopen(KEGG_FTP_PATH + filename), 
                _join(filename))
    
class KEGGGenomeEntry(KEGGDBEntry):
    FIELDS = ["ENTRY", "NAME", "DEFINITION", "ANNOTATION", "TAXONOMY", "DATA_SOURCE",
              "ORIGINAL_DB", "CHROMOSOME", "STATISTICS", "REFERENCE"]
    MULTIPLE_FIELDS = ["REFERENCE"]
    
    def org_code(self):
        if self.name is not None:
            return self.name.split(",")[0]
        else:
            return self.entry.split()[0]
        
    def entry_key(self):
        return self.org_code()
    
entry_decorate(KEGGGenomeEntry)
        
class KEGGGenome(KEGGDataBase):
    Entry = KEGGGenomeEntry
    VERSION = "v1.1"
    TAXID_MAP = {"4932" : "559292"}
    def name(self):
        return super(KEGGGeneome, self).name()
    
    @loads("KEGG", "kegg_genome.tar.gz", version=VERSION)
    def load(self, file=None):
        filename = _join("genes", "genome")
        super(KEGGGenome, self).load(filename)
    
    @classmethod
    def common_organisms(cls):
#        genome = KEGGGenome()
#        id_map = {"562":"511145", "2104":"272634", "5833":"36329", "4896":"284812", "11103":None, "4754":None, "4577":None}
#        return [genome.get(id_map.get(taxid, taxid)).entry_key() for taxid in obiTaxonomy.common_taxids() if id_map.get(taxid, taxid) is not None]
        return ['ath', 'bta', 'cel', 'cre', 'dre', 'ddi', 'dme', 'eco', 'hsa', 'mmu', 'mpn', 'osa',
                'pfa', 'rno', 'sce', 'spo', 'zma', 'xla']
        
    @classmethod
    def essential_organisms(cls):
        genome = KEGGGenome()
        
        return [genome.get(cls.TAXID_MAP.get(taxid, taxid)).entry_key() for taxid in obiTaxonomy.essential_taxids()]
            
    def get(self, key, *args):
        if key in self.TAXID_MAP:
            key = self.TAXID_MAP[key]
        if key not in self:
            keys = self.search(key)
            keys = [k for k in keys if key in self[k].name or key in self[k].taxonomy]
            key = keys.pop(0) if keys else key
        return super(KEGGGenome, self).get(key, *args)
    
    def search(self, string, relevance=False):
        if string in self.TAXID_MAP:
            string = self.TAXID_MAP[string]
            
        def match(entry, string):
            rel = 0
            if string in entry.entrytext:
                weight_f = [(entry.definition or "", 2), (entry.name or "", 4), (entry.taxonomy or "", 8), (entry.entry_key(), 16)]
                rel += 1 + sum([w for text, w in weight_f if string in text]) + (16 if entry.entry_key() in self.common_organisms() else 0)
            return rel
        matched = sorted(zip([match(entry, string) for entry in self.entrys], self.entrys), reverse=True)
        return [(entry.entry_key(), rel) if relevance else entry.entry_key() for rel, entry in matched if rel]
            
    @classmethod
    @downloader
    def download(cls, file=None):
        return (urllib2.urlopen(KEGG_FTP_PATH + "/genes/genome"), 
                _join("genes", "genome") if file is None else file)
        
borg_class(KEGGGenome)
        
class KEGGCompoundEntry(KEGGDBEntry):
    FIELDS = ["ENTRY", "NAME", "FORMULA", "MASS", "REMARK", "REACTION", "PATHWAY",
              "ENZYME", "DBLINKS", "ATOM", "BOND"]
    def entry_key(self):
        return "cpd:" + super(KEGGCompoundEntry, self).entry_key()
    
entry_decorate(KEGGCompoundEntry)
    
class KEGGCompounds(KEGGDataBase):
    Entry = KEGGCompoundEntry
    FILENAME = "%(path)s/ligand/compound/compound"
    
    @loads("KEGG", "kegg_ligand.tar.gz")
    def load(self, file=None):
        super(KEGGCompounds, self).load(file)
        
    @classmethod
    @downloader
    def download(cls, file=None):
        return (urllib2.urlopen(KEGG_FTP_PATH + "/ligand/compound/compound"),
                _join("ligand", "compound", "compound") if file is None else file)
    
borg_class(KEGGCompounds)
    
class KEGGEnzymeEntry(KEGGDBEntry):
    FIELDS = ["ENTRY", "NAME", "CLASS", "SYSNAME", "REACTION", "ALL_REAC", "SUBSTRATE",
              "PRODUCT", "COFACTOR", "COMMENT", "REFERENCE", "PATHWAY", "ORTHOLOGY", 
              "GENES", "STRUCTURES", "DBLINKS"]
    MULTIPLE_FIELDS = ["REFERENCE"]
    def entry_key(self):
        return "EC:" + self.entry.split()[1]
    
entry_decorate(KEGGEnzymeEntry)

class KEGGEnzymes(KEGGDataBase):
    Entry = KEGGEnzymeEntry
    FILENAME = "%(path)s/ligand/enzyme/enzyme"

    @loads("KEGG", "kegg_ligand.tar.gz")
    def load(self, file=None):
        super(KEGGEnzymes, self).load(file)
    
    @classmethod
    @downloader
    def download(cls, file=None):
        return (urllib2.urlopen(KEGG_FTP_PATH + "ligand/enzyme/enzyme"),
                _join("ligand", "enzyme", "enzyme"))
    
borg_class(KEGGEnzymes)

class KEGGReactionEntry(KEGGDBEntry):
    FIELDS = ["ENTRY", "NAME", "DEFINITION", "EQUATION", "COMMENT", "RPAIR", "PATHWAY", 
              "ENZYME", "ORTHOLOGY"]
    def entry_key(self):
        return "rn:" + super(KEGGReactionEntry, self).entry_key()
    
entry_decorate(KEGGReactionEntry)

class KEGGReactions(KEGGDataBase):
    Entry = KEGGReactionEntry
    FILENAME = "%(path)s/ligand/reaction/reaction"
    
    @loads("KEGG", "kegg_ligand.tar.gz")
    def load(self, file=None):
        super(KEGGReactions, self).load(file)
        
    @classmethod
    @downloader
    def download(cls, file=None):
        return (urllib2.urlopen(KEGG_FTP_PATH + "ligand/reaction/reaction"),
                _join("ligand", "reaction", "reaction"))
        
borg_class(KEGGReactions)

class KEGGKOEntry(KEGGDBEntry):
    FIELDS = ["ENTRY", "NAME", "DEFINITION", "CLASS", "DBLINKS", "GENES"]
    
    def entry_key(self):
        return "ko:" + super(KEGGKOEntry, self).entry_key()
    
entry_decorate(KEGGKOEntry)
    
class KEGGKO(KEGGDataBase):
    Entry = KEGGKOEntry
    FILENAME = "%(path)s/genes/ko"
    
    @loads("KEGG", "kegg_orthology.tar.gz")
    def load(self, file=None):
        super(KEGGKO, self).load(file)
        
    @classmethod
    @downloader
    def download(cls, file=None):
        return (urllib2.urlopen(KEGG_FTP_PATH + "genes/ko"),
                _join("genes", "ko"))
        
borg_class(KEGGKO)

class KEGGBriteEntry(object):
    _search_re = {"ids": re.compile('(?P<ids>\[.*:.*\])'),
                  "title": re.compile(r'(<[Bb]>)?(?P<title>\b[a-zA-Z0-9_/\s,;:.+=\-\[\]{}\(\)]+?)(?(1)</[Bb]>)$'),
                  "links": re.compile('(?P<links><a href=".+?">.*?</a>)')}
    def __init__(self, line, entrys=None):
        self.entrys = [] #entrys if entrys is not None else []
        self.line = line[1:].strip()
        for name, re in self._search_re.items():
            search = re.search(self.line)
            setattr(self, name, search.group(name) if search else None)
#        self.identifiers = self.groups.get("ids")
#        self.title = self.groups. # TODO: parse line (html, kegg identifiers ...)

    def __iter__(self):
        return iter(self.entrys)

class KEGGBrite(KEGGBriteEntry):
    VERSION = "v1.0"
    def __init__(self, brite_id):
        super(KEGGBrite, self).__init__("")
        self.load(brite_id)
        
    @loads("KEGG", "kegg_brite.tar.gz")
    def load(self, brite_id):
        file = self.filename(brite_id)
        lines = open(file, "rb").read().split("\n!\n")[1].splitlines()
        def collect(lines, depth, collection):
            while lines:
                line = lines[0]
                if line.startswith("#"):
                    lines.pop(0)
                elif line.startswith(depth) and len(line.strip()) > 1:
                    collection.append(KEGGBriteEntry(lines.pop(0))) 
                elif line[0] > depth:
                    collect(lines, line[0], collection[-1].entrys)
                elif line[0] < depth:
                    return
                else:
                    lines.pop(0)
                        
        collect([line for line in lines if not line.startswith("#") and len(line) > 1], "A", self.entrys)
    
    @classmethod
    @defaultpath
    def filename(cls, brite_id, path=DEFAULT_DATABASE_PATH):
        if path is None:
            path = "%(path)s"
        if brite_id.startswith("br"):
            return ("%(path)s/brite/br/" + brite_id + ".keg") % dict(path=path)
        elif brite_id.startswith("ko"):
            return ("%(path)s/brite/ko/" + brite_id + ".keg") % dict(path=path)
        else:
            org = brite[:-5]
            return ("%(path)s/brite/organisms/" + org + "/" + brite_id + ".keg") % dict(path=path)
        
    @classmethod
    @downloader
    def download(cls, brite_id):
        filename = cls.filename(brite_id, path="")
        return (urllib2.urlopen(cls.filename(brite_id, path=KEGG_FTP_PATH.rstrip("/"))),
                cls.filename(brite_id))

class KEGGOrganism(object):
    VERSION = "v2.0"
    DOMAIN = "KEGG"
    def __init__(self, org, genematcher=None):
        self.load(org)
        self.genematcher = genematcher
        
#    @loads("KEGG", lambda self, org: "kegg_genes_%s.tar.gz" % self.organism_name_search(org))
    def load(self, org):
        org = self.organism_name_search(org)
        self.org_code = org
    
    @property
#    @deprecated("Use org_code instead")
    def org(self):
        return self.org_code
    
    @cached_property
    def genes(self):
        return KEGGGenes(self.org_code)
        
    @cached_property
    def gene_aliases(self):
        return self.genes.gene_aliases()
    
    def pathways(self, with_ids=None):
        return [pathway for pathway, values in KEGGPathway.list(self.org_code).items() if all([id in values for id in (with_ids or [])])]
    
#    @deprecated("Use KEGGOrganism.pathways instead")
    def list_pathways(self):
        return self.pathways()
    
#    @deprecated
    def get_linked_pathways(self, pathway_id):
        return KEGGPathway(pathway_id).genes()
    
    def enzymes(self, genes=None):
        enzymes = KEGGEnzymes()
        return [enzyme.entry_key() for enzyme in enzymes.itervalues() if enzyme.genes and self.org_code in enzyme.genes]
    
    def get_enriched_pathways(self, genes, reference=None, prob=obiProb.Binomial(), callback=None):
        """Return a dictionary with enriched pathways ids as keys and (list_of_genes, p_value, num_of_reference_genes) tuples as items."""
        allPathways = defaultdict(lambda :[[], 1.0, []])
        from Orange.orng import orngMisc
        milestones = orngMisc.progressBarMilestones(len(genes), 100)
        for i, gene in enumerate(genes):
            pathways = self.pathways([gene])
            for pathway in pathways:
                allPathways[pathway][0].append(gene)
            if callback and i in milestones:
                callback(i*100.0/len(genes))
        reference = set(reference if reference is not None else self.genes.keys())
        for p_id, entry in allPathways.items():
            entry[2].extend(reference.intersection(KEGGPathway(p_id).genes()))
            entry[1] = prob.p_value(len(entry[0]), len(reference), len(entry[2]), len(genes))
        return dict([(pid, (genes, p, len(ref))) for pid, (genes, p, ref) in allPathways.items()])
    
    def get_genes_by_enzyme(self, enzyme):
        enzyme = KEGGEnzymes()[enzyme]
        return enzyme.genes.get(self.org_code, []) if enzyme.genes else []
    
    def get_genes_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).genes()
    
    def get_enzymes_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).enzymes()
    
    def get_compounds_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).compounds()
    
    def get_pathways_by_genes(self, gene_ids):
        gene_ids = set(gene_ids)
        pathways = [self.genes[id].pathway for id in gene_ids if self.genes[id].pathway]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways if gene_ids.issubset(KEGGPathway(id).genes())] 
    
    def get_pathways_by_enzymes(self, enzyme_ids):
        enzyme_ids = set(enzyme_ids)
        pathways = [KEGGEnzymes()[id].pathway for id in enzyme_ids]
        pathwats = reduce(set.union, pathways, set())
        return [id for id in pathways if enzyme_ids.issubset(KEGGPathway(id).enzymes())]
    
    def get_pathways_by_compounds(self, compound_ids):
        compound_ids = set(compound_ids)
        pathways = [KEGGCompounds()[id].pathway for id in compound_ids]
        pathwats = reduce(set.union, pathways, set())
        return [id for id in pathways if compound_ids.issubset(KEGGPathway(id).compounds())]
    
    def get_enzymes_by_compound(self, compound_id):
        return KEGGCompound()[compound_id].enzyme
    
    def get_enzymes_by_gene(self, gene_id):
        return self.genes[gene_id].enzymes
    
    def get_compounds_by_enzyme(self, enzyme_id):
        return self._enzymes_to_compounds.get(enzyme_id)
    
    def get_unique_gene_ids(self, genes, caseSensitive=True):
        """Return a tuple with three elements. The first is a dictionary mapping from unique gene
        ids to gene names in genes, the second is a list of conflicting gene names and the third is a list
        of unknown genes.
        """
        unique, conflicting, unknown = {}, [], []
        for gene in genes:
            names = self.genematcher.match(gene)
            if len(names) == 1:
                unique[names[0]] = gene
            elif len(names) == 0:
                unknown.append(gene)
            else:
                conflicting.append(gene)
        return unique, conflicting, unknown
    
    @property
    @cached_method
    def _enzymes_to_compounds(self):
        dd = {}
        for val in KEGGCompounds().values():
            dd.update(dict.fromkeys(val.enzymes, val.entry_key()))
        return dd
    
    def _set_genematcher(self, genematcher):
        setattr(self, "_genematcher", genematcher)
        
    def _get_genematcher(self):
        if getattr(self, "_genematcher", None) == None:
            from . import obiGene
            if self.org_code == "ddi":
                self._genematcher = obiGene.matcher([obiGene.GMKEGG(self.org_code), obiGene.GMDicty(),
                                                     [obiGene.GMKEGG(self.org_code), obiGene.GMDicty()]])
            else:
                self._genematcher = obiGene.matcher([obiGene.GMKEGG(self.org_code)])
            self._genematcher.set_targets(self.genes.keys())
        return self._genematcher
    genematcher = property(_get_genematcher, _set_genematcher)
    
    def get_genes(self):
        return self.genes
    
    @classmethod
    def organism_name_search(cls, name):
        genome = KEGGGenome()
        if name not in genome:
            ids = genome.search(name)
            if not ids:
                from . import obiTaxonomy
                ids = obiTaxonomy.search(name)
                ids = [id for id in ids if genome.search(id)]
            name = ids.pop(0) if ids else name
            
        try:
            return KEGGGenome()[name].entry_key()
        except KeyError:
            raise ValueError("Organims with name='%s' not found in KEGG." % name)
        
    @classmethod
    def organism_version(cls, name):
        name = cls.organism_name_search(name)
        orngServerFiles.localpath_download("KEGG", "kegg_genes_%s.tar.gz" % name)
        return orngServerFiles.info("KEGG", "kegg_genes_%s.tar.gz" % name)["datetime"]
              
class KEGGPathway(object):
    PNG_FILENAME = "%(path)s/pathway/organisms/%(org)s/%(org)s%(map_id)s.png"
    KGML_FILENAME = "%(path)s/xml/kgml/metabolic/organisms/%(org)s/%(org)s%(map_id)s.xml"
    VERSION = "v2.0"
    
    class entry(object):
        def __init__(self, dom_element):
            self.__dict__.update(dom_element.attributes.items())
            self.graphics = ()
            self.components = []
            self.graphics = dict(dom_element.getElementsByTagName("graphics")[0].attributes.items())
            self.components = [node.getAttribute("id") for node in dom_element.getElementsByTagName("component")]
    class reaction(object):
        def __init__(self, dom_element):
            self.__dict__.update(dom_element.attributes.items())
            self.substrates = [node.getAttribute("name") for node in dom_element.getElementsByTagName("substrate")]
            self.products = [node.getAttribute("name") for node in dom_element.getElementsByTagName("product")]
            
    class relation(object):
        def __init__(self, dom_element):
            self.__dict__.update(dom_element.attributes.items())
            self.subtypes = [node.attributes.items() for node in dom_element.getElementsByTagName("subtype")]
        
    def __init__(self, file):
        if not os.path.exists(file):
            path, org, map_id = self.split_pathway_id(file)
#            file = self.KGML_FILENAME % dict(path=DEFAULT_DATABASE_PATH, org=org, map_id=id)
            file = self.filename_kgml(org, map_id)
        if not os.path.exists(file) and os.path.exists(file.replace("metabolic", "non-metabolic")):
            file = file.replace("metabolic", "non-metabolic")
        self.filename = file
#        self.load(file)
        
    @persistent_cached_method
    def pathway_attributes(self):
        return dict(self.pathway_dom().attributes.items())
    
    @property
    def name(self):
        return self.pathway_attributes().get("name")
    
    @property
    def org(self):
        return self.pathway_attributes().get("org")
    
    @property
    def number(self):
        return self.pathway_attributes().get("number")
    
    @property    
    def title(self):
        return self.pathway_attributes().get("title")
    
    @property
    def image(self):
        return self.pathway_attributes().get("image")
    
    @property
    def link(self):
        return self.pathway_attributes().get("link")
    
    @cached_method
    @loads("KEGG", lambda self: "kegg_pathways_%s.tar.gz" % self.filename.split("/")[-1][:-9])
    def pathway_dom(self):
        return minidom.parse(self.filename).getElementsByTagName("pathway")[0]
    
    @cached_method
    def entrys(self):
        return [self.entry(e) for e in self.pathway_dom().getElementsByTagName("entry")]
    
    @cached_method
    def reactions(self):
        return [self.reaction(e) for e in self.pathway_dom().getElementsByTagName("reaction")]
    
    @cached_method
    def relations(self):
        return [self.relation(e) for e in self.pathway_dom().getElementsByTagName("relation")]
    
    @classmethod
    def split_pathway_id(cls, id):
        path, id = id.split(":") if ":" in id else ("path", id)
        org, id = id[:-5], id[-5:]
        return path, org, id  
        
    def __iter__(self):
        """ Iterate over all elements in the pathway
        """
        return iter(self.all_elements())
    
    def __contains__(self, element):
        """ Retrurn true if element in the pathway
        """
        return element in self.all_elements()
    
    def __getitem__(self, key):
        return
    
    @cached_method
    def all_elements(self):
        """ Return all elements
        """
        return reduce(list.__add__, [self.genes(), self.compounds(), self.enzmes(), self.reactions()], [])
    
    def _get_entrys_by_type(self, type):
        return reduce(set.union, [entry.name.split() for entry in self.entrys() if entry.type == type], set())
    
    @persistent_cached_method
    def genes(self):
        """ Return all genes on the pathway
        """
        return self._get_entrys_by_type("gene")
    
    @persistent_cached_method
    def compounds(self):
        """ Return all compounds on the pathway
        """
        return self._get_entrys_by_type("compound")
    
    @persistent_cached_method
    def enzymes(self):
        """ Return all enzymes on the pathway
        """
        return self._get_entrys_by_type("enzyme")
    
    @persistent_cached_method
    def orthologs(self):
        """ Return all orthologs on the pathway
        """
        return self._get_entrys_by_type("ortholog")
    
    @persistent_cached_method
    def maps(self):
        """ Return all linked maps on the pathway
        """
        return self._get_entrys_by_type("map")
    
    @persistent_cached_method
    def groups(self):
        """ Return all groups on the pathway
        """
        return self._get_entrys_by_type("ortholog")

    def get_image(self):
        """ Return an image of the pathway
        """
        return self.filename_png(self.org, self.number) % dict(path=DEFAULT_DATABASE_PATH)
    
    @persistent_cached_method
    def get_bounding_box_dict(self):
        return dict([(element.id, element.graphics) for element in self.entrys() if element.graphics])
    
    @persistent_cached_method
    def graphics(self, item):
        return [entry.graphics for entry in self.entrys() if item in entry.name and entry.graphics]
        
    @classmethod
    @partial(cached_method, cache_name="_cls_cached_method_cache_")
    @loads("KEGG", lambda cls, org: "kegg_pathways_%s.tar.gz" % org)
    def list(cls, org):
        file = cls.directory_png(org) + org + ".list"
        data = [line.split() for line in open(file, "rb").read().splitlines()]
        return reduce(lambda dict, line: dict[line[0]].update(line[1:]) or dict, data, defaultdict(set))
    
    @classmethod
    @defaultpath
    def filename_kgml(cls, org, map_id, path=DEFAULT_DATABASE_PATH):
        path = "%(path)s" if path is None else path
        if org in ["map", "ec"]:
            return "%(path)s/xml/kgml/metabolic/ec/ec%(map_id)s.xml" % dict(org=org, map_id=map_id, path=path)
        elif org == "ko":
            return "%(path)s/xml/kgml/metabolic/ko/ko%(map_id)s.xml" % dict(org=org, map_id=map_id, path=path)
        else:
            return "%(path)s/xml/kgml/metabolic/organisms/%(org)s/%(org)s%(map_id)s.xml" % dict(org=org, map_id=map_id, path=path)
        
    @classmethod
    @defaultpath
    def filename_png(cls, org, map_id, path=DEFAULT_DATABASE_PATH):
        path = "%(path)s" if path is None else path 
        if org in ["map", "ec", "ko", "rn"]:
            return "%(path)s/pathway/%(org)s/%(org)s%(map_id)s.png" % dict(org=org, map_id=map_id, path=path)
        else:
            return "%(path)s/pathway/organisms/%(org)s/%(org)s%(map_id)s.png" % dict(org=org, map_id=map_id, path=path)
        
    @classmethod
    @defaultpath
    def directory_png(cls, org, path=DEFAULT_DATABASE_PATH):
        return cls.filename_png(org, "", path=path).rsplit("/", 1)[0] + "/" 
    
    @classmethod
    @defaultpath
    def directory_kgml(cls, org, path=DEFAULT_DATABASE_PATH):
        return cls.filename_kgml(org, "", path=path).rsplit("/", 1)[0] + "/"
        
    @classmethod
    def pathways(cls, org):
        file = cls.directory_png(org) + org + ".list"
        pathways = [line.split()[0] for line in open(file, "rb").read().splitlines() if line.strip()]
        return sorted(set(pathways))
        
    @classmethod
    def download(cls, pathway, target=None):
        """ Download the pathway (xml and png files) to target directory
        """
        org = pathway[:3]
        import urllib2
        xml_data = urllib2.urlopen("ftp://ftp.genome.jp/pub/kegg/xml/kgml/metabolic/organisms/%s/%s.xml" % (org, pathway)).read()
        open(os.path.join(target, pathway + ".xml"), "wb").write(xml_data)
        png_data = urllib2.urlopen("ftp://ftp.genome.jp/pub/kegg/pathway/organisms/%s/%s.png" % (org, pathway)).read()
        open(os.path.join(target, pathway + ".png"), "wb").write(png_data)
        
    @classmethod
    @downloader
    def download_pathways(cls, org):
        png_path = cls.directory_png(org, path=None).replace("%(path)s/", "")
        xml_path = cls.directory_kgml(org, path=None).replace("%(path)s/", "")
        data = urllib2.urlopen("ftp://ftp.genome.jp/pub/kegg/" + png_path + org + ".list").read()
        pathways = sorted(set([line.split()[0][5:] for line in data.splitlines() if line]))
        from .obiData import FtpDownloader
        ftp = FtpDownloader("ftp.genome.jp", _join(), "pub/kegg/", numOfThreads=15)
        ftp.massRetrieve([png_path + pathway + ".png" for pathway in pathways])
            
        if org != "map" :
            pathways = ftp.listdir(xml_path)
            pathways = [p for p in pathways if p.endswith(".xml")]
            ftp.massRetrieve([xml_path + pathway for pathway in pathways])
            if org not in ["ec", "rn"]: # Add non-metabolic pathways to the download list
                xml_path = xml_path.replace("metabolic", "non-metabolic")
                pathways = ftp.listdir(xml_path)
                pathways = [p for p in pathways if p.endswith(".xml")] 
                ftp.massRetrieve([xml_path + pathway for pathway in pathways])
        ftp.massRetrieve([png_path + org + ".list"])
        return []
    
persistent_cached_class(KEGGPathway)
        
def organism_name_search(name):
    return KEGGOrganism.organism_name_search(name)

def pathways(org):
    return KEGGPathway.list(org)

def organisms():
    return KEGGOrganism.organisms()

def to_taxid(name):
    genome = KEGGGenome()
    names = genome.search(name)
    if genome[names[0]].taxonomy:
        return re.findall(r"TAX:\s*(\d+)", genome[names[0]].taxonomy)[0]
    else:
        return None #Should this raise an error?

def from_taxid(taxid):
    return KEGGGenome().search(taxid)[0]
    
def test():
    p = KEGGPathway("sce00010.xml")
    print p.genes
    print p.reactions
    print p.compounds
    print p.image
    g = KEGGGenome()
    org = KEGGOrganism("Homo sapiens")
    print list(org.genes)[:10]
    org.gene_aliases
    print org.pathways(with_ids=org.genes.keys()[:5])
    print org.enzymes()
    print org.enriched_pathways(org.genes.keys()[:10])
    print org.genematcher

if __name__=="__main__":
    test()
#    org1 = KEGGOrganism("ddi")
#    org2 = KEGGOrganism("ddi")
#    org2.api = KEGGInterface()
#    tests = [("get_genes", ()),
#             ("get_genes_by_enzyme", ("ec:1.1.1.1",)),
#             ("get_genes_by_pathway", ("path:ddi00010",)),
#             ("get_pathways_by_genes", (["ddi:DDB_0191256"],)),
#             ("get_pathways_by_enzymes", (["ec:1.1.1.1"],)),
#             ("get_pathways_by_compounds", (["cpd:C00001"],)),
#             ("get_linked_pathways", ("path:ddi00010",)),
#             ("list_pathways", ()),
#             ("get_compounds_by_enzyme", ("ec:1.1.1.1",)),
#             ("get_compounds_by_pathway", ("path:ddi00010",)),
#             ("get_enzymes_by_compound", ("cpd:C00001",)),
#             ("get_enzymes_by_pathway", ("path:ddi00010",)),
#             ("get_enzymes_by_gene", ("ddi:DDB_0191256",))]
#    for name, args in tests:
#        s1 = set(getattr(org1, name)(*args))
#        s2 = set(getattr(org2, name)(*args))
#        if s1 and s2:
#            print name
#            print s1-s2
#            print s2-s1
#        else:
#            print name
#            print "both empty"
