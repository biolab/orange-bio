"""
DBGET database
"""
from __future__ import absolute_import

import re

from . import entry
from .entry import fields
from . import api

def iter_take(source_iter, n):
    source_iter = iter(source_iter)
    return [item for _, item in zip(range(n), source_iter)]

def batch_iter(source_iter, n):
    source_iter = iter(source_iter)
    while True:
        batch = iter_take(source_iter, n)
        if batch:
            yield batch
        else:
            break
        
def chain_iter(chains_iter):
    for iter in chains_iter:
        for element in iter:
            yield element

class DBDataBase(object):
    """ A wrapper for DBGET database.
    """
    # ENTRY_TYPE constructor (type)
    ENTRY_TYPE = entry.DBEntry
    
    # Needs to be set in a subclass or object instance 
    DB = None
    
    def __init__(self, **kwargs):
        if not self.DB:
            raise TypeError("Cannot make an instance of abstract base class %r." \
                            % type(self).__name__)
            
        self.api = api.CachedKeggApi()
        self.info = self.api.binfo(self.DB)
        release = self.info.release
        self.api.set_default_release(release)
        self._keys = []
        
    def keys(self):
        return list(self._keys)
    
    def iterkeys(self):
        return iter(self._keys)
    
    def items(self):
        return list(zip(self.keys(), self.batch_get(self.keys())))
    
    def iteritems(self):
        batch_size = 100
        iterkeys = self.iterkeys()
        return chain_iter(zip(batch, self.batch_get(batch))
                          for batch in batch_iter(iterkeys, batch_size))
        
#        return ((key, self.__getitem__(key)) for key in self.iterkeys())
    
    def values(self):
        return self.batch_get(self.keys())
    
    def itervalues(self):
        batch_size = 100
        iterkeys = self.iterkeys()
        return chain_iter(self.batch_get(batch)
                          for batch in batch_iter(iterkeys, batch_size))
        
#        return (self.__getitem__(key) for key in self.iterkeys())
    
    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return default
        
    def has_key(self, key):
        return self.__contains__(key)
    
    def __getitem__(self, key):
        e = self.get_entry(key)
        if e is None:
            raise KeyError(key)
        else:
            return e
    
    def __contains__(self, key):
        return key in set(self.keys())
    
    def __len__(self):
        return len(self.keys())
    
    def __iter__(self):
        return iter(self.keys())
    
    def get_text(self, key):
        key = self._add_db(key)
        return self.api.bget([key])
    
    def get_entry(self, key):
        text = self.get_text(key)
        if not text or text == "None":
            return None
        else:
            return self.ENTRY_TYPE(text)
        
    def find(self, name):
        """ Find ``name`` using BFIND. 
        """
        res = self.api.bfind(self.DB, name).splitlines()
        return [r.split(" ", 1)[0] for r in res]    
        
    def pre_cache(self, keys=None, batch_size=100, progress_callback=None):
        """ Retrive all the entries and cache them locally.
        """
        # TODO do this in multiple threads
    
        if not isinstance(self.api, api.CachedKeggApi):
            raise TypeError("Not an an instance of api.CachedKeggApi")
        
        if batch_size > 100 or batch_size < 1:
            raise ValueError("Invalid batch_size")
        
        if keys is None:
            keys = self.keys()
            
        keys = list(keys)
        start = 0
        while start < len(keys):
            batch = keys[start: start + batch_size]
            batch = map(self._add_db, batch)
            
            self.api.bget(batch)
            
            if progress_callback:
                progress_callback(100.0 * start / len(keys))
                
            start += batch_size
            
    def batch_get(self, keys):
        """ Batch retrieve all entries for keys. This can be
        significantly faster then getting each entry separately
        especially if entries are not yet cached.
        
        """
        entries = []
        batch_size = 100
        keys = list(keys)
        start = 0
        while start < len(keys):
            batch = keys[start: start + batch_size]
            batch = map(self._add_db, batch)
            batch_entries = self.api.bget(batch)
            if batch_entries is not None:
                batch_entries = batch_entries.split("///\n")
                # Remove possible empty last line  
                batch_entries = [e for e in batch_entries if e.strip()]
                entries.extend(map(self.ENTRY_TYPE, batch_entries))
            start += batch_size
            
        return entries
            
    def _add_db(self, key):
        """ Prefix the key with '%(DB)s:' string if not already
        prefixed. 
        """
        if not key.startswith(self.DB + ":"):
            return self.DB + ":" + key
        else:
            return key
        
    @property
    def entries(self):
        return self.values()
    
@entry.entry_decorate
class GenomeEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("DEFINITION", fields.DBDefinitionField),
              ("ANNOTATION", fields.DBSimpleField),
              ("TAXONOMY", fields.DBTaxonomyField),
              ("DATA_SOURCE", fields.DBSimpleField),
              ("ORIGINAL_DB", fields.DBSimpleField),
              ("KEYWORDS", fields.DBSimpleField),
              ("DISEASE", fields.DBSimpleField),
              ("COMMENT", fields.DBSimpleField),
              ("CHROMOSOME", fields.DBFieldWithSubsections),
              ("STATISTICS", fields.DBSimpleField),
              ("REFERENCE", fields.DBReference)]
    
    MULTIPLE_FIELDS = ["REFERENCE"]
    
    def __init__(self, text):
        entry.DBEntry.__init__(self, text)
        
    @property
    def entry_key(self):
        """ Primary entry key used for querying.
        
        .. note:: Unlike most of the other entry types this is the
            first listed 'NAME'.
            
        """
        
        return self.name.split(",", 1)[0]

    @property
    def taxid(self):
        return self.TAXONOMY.taxid
            
    def org_code(self):
        if self.name is not None:
            return self.name.split(",")[0]
        else:
            return self.entry.split(" ")[0]
        

class Genome(DBDataBase):
    DB = "genome"
    ENTRY_TYPE = GenomeEntry
    
    # For obiTaxonomy.common_taxids mapping
    TAXID_MAP = {"562": "511145",   # Escherichia coli K-12 MG1655
                 "2104": "272634",  # Mycoplasma pneumoniae M129 
                 "4530": "39947",   # Oryza sativa ssp. japonica cultivar Nipponbare (Japanese rice)
                 "4932" : "559292", # Saccharomyces cerevisiae S288C
                 "4896": "284812",  # Schizosaccharomyces pombe 972h-
                 }
    
    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [org.entry_id for org in self.api.list_organisms()]
    
    def _key_to_gn_entry_id(self, key):
        res = self.find(key)
        if len(res) == 0:
            raise KeyError("Unknown key")
        elif len(res) > 1:
            raise ValueError("Not a unique key")
        else:
            return res[0]
    
    @classmethod
    def common_organisms(cls):
        return ['ath', 'bta', 'cel', 'cre', 'dre', 'ddi',
                'dme', 'eco', 'hsa', 'mmu', 'mpn', 'osa',
                'pfa', 'rno', 'sce', 'spo', 'zma', 'xla']
        
    @classmethod
    def essential_organisms(cls):
        return ['ddi', 'dme', 'hsa', 'mmu', 'sce']
    
    def search(self, string, relevance=False):
        """ Search the genome database for string using ``bfind``.
        """
        if relevance:
            raise NotImplementedError("relevance is no longer supported")
        if string in self.TAXID_MAP:
            string = self.TAXID_MAP[string]
            
        res = self.api.bfind(self.DB, string)
        if not res:
            return []
        
        res = res.splitlines()
        res = [r.split(",", 1)[0] for r in res]
        res = [r.split(" ", 1)[1] for r in res]
        return res
    
    
@entry.entry_decorate
class GeneEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("DEFINITION", fields.DBDefinitionField),
              ("ORGANISM", fields.DBSimpleField),
              ("ORTHOLOGY", fields.DBSimpleField),
              ("DRUG_TARGET", fields.DBSimpleField),
              ("PATHWAY", fields.DBPathway),
              ("MODULE", fields.DBSimpleField),
              ("DISEASE", fields.DBSimpleField),
              ("CLASS", fields.DBSimpleField),
              ("POSITION", fields.DBSimpleField),
              ("MOTIF", fields.DBSimpleField),
              ("DBLINKS", fields.DBDBLinks),
              ("STRUCTURE", fields.DBSimpleField),
              ("AASEQ", fields.DBAASeq),
              ("NTSEQ", fields.DBNTSeq)]
    
    def aliases(self):
        return [self.entry_key] + (self.name.split(",") if self.name else []) + [link[1][0] for link in self.dblinks.items() if self.dblinks]

    @property
    def alt_names(self):
        """ For backwards compatibility.
        """
        return self.aliases()
  
class Genes(DBDataBase):
    DB = None # Needs to be set in __init__ 
    ENTRY_TYPE = GeneEntry
    
    def __init__(self, org_code):
        self.DB = org_code
        self.org_code = org_code
        DBDataBase.__init__(self)
        self._keys = self.api.get_genes_by_organism(org_code)
        
    def gene_aliases(self):
        aliases = {}
        for entry in self.itervalues():
            aliases.update(dict.fromkeys(entry.aliases(), self.org_code + ":" + entry.entry_key()))
        return aliases
    

@entry.entry_decorate
class CompoundEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("FORMULA", fields.DBSimpleField),
              ("MASS", fields.DBSimpleField),
              ("REMARK", fields.DBSimpleField),
              ("REACTION", fields.DBSimpleField),
              ("PATHWAY", fields.DBPathway),
              ("ENZYME", fields.DBSimpleField),
              ("DBLINKS", fields.DBDBLinks),
              ("ATOM", fields.DBSimpleField),
              ("BOND", fields.DBSimpleField)
              ]
    
    
class Compounds(DBDataBase):
    DB = "cpd"
    ENTRY_TYPE = CompoundEntry
    
    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [] # All keys are not available


@entry.entry_decorate    
class ReactionEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("DEFINITION", fields.DBDefinitionField),
              ("EQUATION", fields.DBSimpleField),
              ("ENZYME", fields.DBSimpleField)
              ]
    
class Reactions(DBDataBase):
    DB = "rn"
    ENTRY_TYPE = ReactionEntry
    
    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [] # All keys are not available
         
class Brite(DBDataBase):
    DB = "br"
    
class Disease(DBDataBase):
    DB = "ds"
        
class Drug(DBDataBase):
    DB = "dr"
    
@entry.entry_decorate
class EnzymeEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("CLASS", fields.DBSimpleField),
              ("SYSNAME", fields.DBSimpleField),
              ("REACTION", fields.DBSimpleField),
              ("ALL_REAC", fields.DBSimpleField),
              ("SUBSTRATE", fields.DBSimpleField),
              ("PRODUCT", fields.DBSimpleField),
              ("COMMENT", fields.DBSimpleField),
              ("REFERENCE", fields.DBReference),
              ("PATHWAY", fields.DBPathway),
              ("ORTHOLOGY", fields.DBSimpleField),
              ("GENES", fields.DBSimpleField),
              ("DBLINKS", fields.DBDBLinks)
              ]
    
    MULTIPLE_FIELDS = ["REFERENCE"]
    
class Enzymes(DBDataBase):
    DB = "ec"
    ENTRY_TYPE = EnzymeEntry
    
    
@entry.entry_decorate
class OrthologyEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("CLASS", fields.DBSimpleField),
              ("DBLINKS", fields.DBDBLinks),
              ("GENES", fields.DBSimpleField),
              ]
    
class Orthology(DBDataBase):
    DB = "ko"
    ENTRY_TYPE = OrthologyEntry
    
    
@entry.entry_decorate
class PathwayEntry(entry.DBEntry):
    FIELDS = [("ENTRY", fields.DBEntryField),
              ("NAME", fields.DBNameField),
              ("DESCRIPTION", fields.DBSimpleField),
              ("CLASS", fields.DBSimpleField),
              ("PATHWAY_MAP", fields.DBPathwayMapField),
              ("DISEASE", fields.DBSimpleField),
              ("DRUG", fields.DBSimpleField),
              ("DBLINKS", fields.DBDBLinks),
              ("ORGANISM", fields.DBSimpleField),
              ("GENE", fields.DBGeneField),
              ("ENZYME", fields.DBEnzymeField),
              ("COMPOUND", fields.DBCompoundField),
              ("REFERENCE", fields.DBReference),
              ("REL_PATHWAY", fields.DBSimpleField),
              ("KO_PATHWAY", fields.DBSimpleField),
              ]
    
    MULTIPLE_FIELDS = ["REFERENCE"]
    
    @property
    def gene(self):
        if hasattr(self, "GENE"):
            genes = self.GENE._convert()
        else:
            return None
        
        org = self.organism
        org_prefix = ""
        if org:
            match = re.findall(r"\[GN:([a-z]+)\]", org)
            if match:
                org_prefix = match[0] + ":"
        genes = [org_prefix + g for g in genes]
        return genes 
    
class Pathways(DBDataBase):
    DB = "path"
    ENTRY_TYPE = PathwayEntry
    
    def __init__(self):
        DBDataBase.__init__(self)
    