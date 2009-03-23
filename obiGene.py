
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
            

    @classmethod    
    def load(cls, file):
        """ A class method that loads gene info from file
        """
        if type(file) in [str, unicode]:
            file = open(file, "rb")
        return cls((line.split("\t", 3)[1], line) for line in file.read().split("\n") if line.strip() and not line.startswith("#"))
        
    def get_info(self, gene_id):
        """ Search and return the GeneInfo object for gene_id
        """
        return self[gene_id]
        
    def __call__(self, name):
        """ Search and return the GeneInfo object for gene_id
        """
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
        import urllib2, gzip
        from cStringIO import StringIO
        data = StringIO(urllib2.urlopen("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz").read())
        info = NCBIGeneInfo.load(gzip.GzipFile(None, "rb", fileobj=data))
        return info
        

if __name__ == "__main__":
    info = NCBIGeneInfo("9606")
    gi = info(list(info)[0])
    print gi.tax_id, gi.synonyms, gi.dbXrefs, gi.symbol_from_nomenclature_authority, gi.full_name_from_nomenclature_authority
    