import orngServerFiles
import sys, os
import urllib2
import re

from collections import defaultdict

class disease(object):
    regex = re.compile(r'(?P<name>.*?),? (?P<id>[0-9]{3,6} )?(?P<m1>\([123?]\) )?(?P<m2>\([123?]\) )? *$')
    __slots__ = ["name", "id", "mapping"]
    def __init__(self, morbidmap_line):
        string = morbidmap_line.split("|", 1)[0]
        match = self.regex.match(string)
#        print string
#        print match.groups()
        self.name, self.id, self.mapping = [s.strip() if s else s for s in match.groups()[:3]]
        if match.group("m2"):
            self.mapping += " " + match.group("m2").strip()
        
class OMIM(object):
    VERSION = 1
    DEFAULT_DATABASE_PATH = orngServerFiles.localpath("OMIM")
    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path else self.DEFAULT_DATABASE_PATH
        self.load()
    
    @classmethod
    def download_from_NCBI(cls, file=None):
        data = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/repository/OMIM/morbidmap").read()
        if file is None:
            if not os.path.exists(cls.DEFAULT_DATABASE_PATH):
                os.mkdir(cls.DEFAULT_DATABASE_PATH)
            file = open(os.path.join(cls.DEFAULT_DATABASE_PATH, "morbidmap"), "wb")
        elif type(file) in [str, unicode]:
            file = open(file, "wb")
        file.write(data)
        
    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            omim = OMIM()
            cls._shared_dict = omim.__dict__
        instance = OMIM.__new__(OMIM)
        instance.__dict__ = cls._shared_dict
        return instance 
    
    def load(self):
        orngServerFiles.localpath_download("OMIM", "morbidmap")
        lines = open(os.path.join(self.local_database_path, "morbidmap")).read().split("\n")
        self._disease_dict = dict([(disease(line), line) for line in lines if line])
        
    def diseases(self):
        return self._disease_dict.keys()
    
    def genes(self):
        return sorted(set(reduce(list.__add__, [self.disease_genes(disease) for disease in self.diseases()], [])))
    
    def disease_genes(self, disease):
        return self._disease_dict[disease].split("|")[1].split(", ")
    
    def gene_diseases(self):
        d = defaultdict(set)
        for disease, genes in [(disease, self.disease_genes(disease)) for disease in self.diseases()]:
            for gene in genes:
                d[gene].add(disease)
        return d
def diseases():
    return OMIM.get_instance().diseases()
        
def genes():
    return OMIM.get_instance().genes()

def disease_genes(disease):
    return OMIM.get_instance().disease_genes(disease)

def genes_diseases():
    return OMIM.get_instance().gene_diseases()
