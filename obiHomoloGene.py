import orngServerFiles
import sys, os
import urllib2

from collections import defaultdict

class _homolog(object):
    __slots__ = ["group_id", "taxonomy_id", "gene_id", "gene_symbol"]
    def __init__(self, homolog_line):
        for attr, val in zip(self.__slots__, homolog_line.split("\t")):
            setattr(self, attr, val)
    
class HomoloGene(object):
    DEFAULT_DATABASE_PATH = orngServerFiles.localpath("HomoloGene")
    VERSION = 1
    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path else self.DEFAULT_DATABASE_PATH
        self.load()

    @classmethod
    def download_from_NCBI(cls, file=None):
        data = urllib2.urlopen("ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current/homologene.data").read()
        if file is None:
            try:
                os.mkdir(orngServerFiles.localpath("HomoloGene"))
            except OSError:
                pass
            file = open(orngServerFiles.localpath("HomoloGene", "homologene.data"), "wb")
        elif type(file) in [str, unicode]:
            file = open(file, "wb")
        file.write(data)
        file.flush()
        
    @classmethod    
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            h = HomoloGene()
            cls._shared_dict = h.__dict__
        h = cls.__new__(cls)
        h.__dict__ = cls._shared_dict
        return h
    
    def load(self):
        orngServerFiles.localpath_download("HomoloGene", "homologene.data")
        lines = open(orngServerFiles.localpath("HomoloGene", "homologene.data"), "rb").read().split("\n")[:-1]
        self._homologs = {}
        self._homologs_by_group = defaultdict(list)
        self._homologs = dict([((h.taxonomy_id, h.gene_symbol), h) for h in [_homolog(line) for line in lines]])
#        for line in lines:
#            h = _homolog(line)
#            self._homologs[h.taxonomy_id, h.gene_symbol] = h
#            self._homologs_by_group[h.group_id].append(h)
        
    def all_genes(self, taxid=None):
        return [homolog.gene_symbol for (tid, id), homolog in self._homologs.iteritems() if tid == taxid]
    
    def homologs(self, gene, taxid):
        group = self._homologs.get((taxid, gene), _homolog("")).group_id
        homologs = [h for h in self._homologs.itervalues() if h.group_id == group] #self._homologs_by_group[group]
        return [(h.taxonomy_id, h.gene_symbol) for h in homologs]
        
    def homolog(self, gene, taxid, homolotaxid):
        homologs = dict(self.homologs(gene, taxid))
        return homologs.get(homolotaxid, None)
        
def all_genes(taxid):
    """ Return a set of all genes for organism with taxid
    """
    return HomoloGene.get_instance().all_genes(taxid)

def homologs(genename, taxid):
    """ Return a list of homologs (taxid, genename) for a homolog group that gene, taxid belong to 
    """ 
    return HomoloGene.get_instance().homologs(genename, taxid)

def homolog(genename, taxid, homolotaxid):
    """ Return a homolog of genename, taxid in organism with holomotaxid or None if homolog does not exist.
    """
    return HomoloGene.get_instance().homolog(genename, taxid, homolotaxid)