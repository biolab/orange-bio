
ncbi_geneinfo_tags = ("tax_id", "gene_id", "symbol", "locus_tag", "synonyms",
                      "dbXrefs", "chromosome", "map_location", "description", "type",
                      "symbol_from_nomenclature_authority", "full_name_from_nomenclature_authority",
                      "nomenclature_status", "other_designations", "modification_date")

ncbi_multiple_cardinality_tags = ("synonyms", "dbXrefs", "other_designations")

import os
import obiTaxonomy
import orngServerFiles
from obiTaxonomy import TextDB
from weakref import WeakValueDictionary

default_database_path = orngServerFiles.localpath("NCBI_geneinfo")

class GeneInfo(object):
    __slots__ = ("tax_id", "gene_id", "symbol", "locus_tag", "synonyms",
                 "dbXrefs", "chromosome", "map_location", "description", "type",
                 "symbol_from_nomenclature_authority", "full_name_from_nomenclature_authority",
                 "nomenclature_status", "other_designations", "modification_date")
    def __init__(self, line):
        line = line.split("\t")
##        line = line[:1] + line[2:]
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
##    _object_cache = WeakValueDictionary()
    _object_cache = {}
    def __init__(self, *args, **kwargs):
        """ An object for accessing NCBI gene info
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
##            self.taxid = args[0]
            self.taxid = taxids.pop()
            if not os.path.exists(orngServerFiles.localpath("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)):
                orngServerFiles.download("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)
            file = open(orngServerFiles.localpath("NCBI_geneinfo", "gene_info.%s.db" % self.taxid), "rb")
            self.update(dict((line.split("\t", 3)[1], line) for line in file.read().split("\n") if line.strip() and not line.startswith("#")))
##            if self.taxid  not in self._object_cache:
##                if not os.path.exists(orngServerFiles.localpath("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)):
##                    orngServerFiles.download("NCBI_geneinfo", "gene_info.%s.db" % self.taxid)
##                self._object_cache[self.taxid] = self.load(orngServerFiles.localpath("NCBI_geneinfo", "gene_info.%s.db" % self.taxid))
##                
##            self.__dict__ = self._object_cache[self.taxid].__dict__
        else:
            dict.__init__(self, *args, **kwargs)
            

    @classmethod    
    def load(cls, file):
        if type(file) in [str, unicode]:
            file = open(file, "rb")
##        self._data = dict([(line.split("\t", 3)[1], line) for line in file.read().split("\n") if line])
##        self.update(dict([(line.split("\t", 3)[1], line) for line in file.read().split("\n") if line]))
        return cls((line.split("\t", 3)[1], line) for line in file.read().split("\n") if line.strip() and not line.startswith("#"))
        
    def get_info(self, id):
##        return GeneInfo(self._data[id])
        return self[id]
        
    def __call__(self, name):
        translate = lambda a:a
        id = translate(name)
        return self.get_info(id)

    def __getitem__(self, key):
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
##        if type(file) in [unicode, str]:
##            file = open(file, "wb")
        import urllib2, gzip
        from cStringIO import StringIO
        data = StringIO(urllib2.urlopen("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz").read())
        info = NCBIGeneInfo.load(gzip.GzipFile(None, "rb", fileobj=data))
        return info
        
def prepare_gene_info(file):
    info = NCBIGeneInfo.load(file)
    taxids = ["9606"]#obiTaxonomy.common_organisms()
    genes = dict([(taxid, []) for taxid in taxids])
    for gi in info.itervalues():
        if gi.tax_id in taxids:
            genes[gi.tax_id].append(gi)

    for taxid, genes in genes.items():
        f = open("gene_info.%s.db" % taxid, "wb")
        f.write("\n".join([str(gi) for gi in sorted(genes, key=lambda gi:int(gi.gene_id))]))
        

if __name__ == "__main__":
##    info = NCBIGeneInfo("9606")
##    gi = info("1")
##    print type(gi), type(info.get("1"))
##    print type(info.values()[0])
    prepare_gene_info("D:/Download/gene_info/Homo_sapiens.gene_info")
##    print gi.tax_id, gi.synonyms, gi.dbXrefs, gi.symbol_from_nomenclature_authority, gi.full_name_from_nomenclature_authority
    