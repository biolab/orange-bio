import sys, os
import urllib2
import shutil
import re
import Orange

from collections import defaultdict

class disease(object):
    """ A disease in the OMIM database.
    """
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
    DEFAULT_DATABASE_PATH = Orange.utils.serverfiles.localpath("OMIM")
    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path is not None else self.DEFAULT_DATABASE_PATH
  
        if self.local_database_path == self.DEFAULT_DATABASE_PATH:
            filename = Orange.utils.serverfiles.localpath_download("OMIM", "morbidmap")
        else:
            filename = os.path.join(self.local_database_path, "morbidmap")

        self.load(filename)
    
    @classmethod
    def download_from_NCBI(cls, file=None):
        if isinstance(file, basestring):
            file = open(file, "wb")
        stream = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/repository/OMIM/ARCHIVE/morbidmap")
        shutil.copyfileobj(stream, file, length=10)
        file.close()

    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            omim = OMIM()
            cls._shared_dict = omim.__dict__
        instance = OMIM.__new__(OMIM)
        instance.__dict__ = cls._shared_dict
        return instance 
    
    def load(self, filename):
        file = open(filename, "rb")
        lines = file.read().splitlines()
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
    """ Return all disease descriptors.
    """
    return OMIM.get_instance().diseases()

def genes():
    """ Return a set of all genes referenced in OMIM.
    """
    return OMIM.get_instance().genes()

def disease_genes(disease):
    """ Return a set of all genes referenced by disease in OMIM.
    """
    return OMIM.get_instance().disease_genes(disease)

def gene_diseases():
    """ Return a dictionary {gene: set(disease_objects for gene), ...}.
    """
    return OMIM.get_instance().gene_diseases()

