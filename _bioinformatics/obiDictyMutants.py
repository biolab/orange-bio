import os
import urllib2
import shutil

from collections import defaultdict

from Orange.orng import orngServerFiles

class DictyMutant(object):

    def __init__(self, mutant_line):
        mutant = mutant_line.split("\t")
        self.id = mutant[0]
        self.descriptor = set(mutant[1].split("/"))
        self.genes = mutant[2].split(" | ")
        self.phenotypes = mutant[3].split(" | ")


class DictyMutants(object):
    VERSION=1
    DEFAULT_DATABASE_PATH = orngServerFiles.localpath("DictyMutants") #use a default folder for storing the genesets

    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path is not None else self.DEFAULT_DATABASE_PATH

        if not os.path.exists(self.local_database_path):
            self.download_from_dictybase(self.local_database_path)

        filename = os.path.join(self.local_database_path, "DictyMutants")
        self.load(filename)

    @classmethod
    def download_from_dictybase(cls, local_database_path=None):
        cls.local_database_path = local_database_path if local_database_path is not None else cls.DEFAULT_DATABASE_PATH

        if not os.path.exists(cls.local_database_path):
            os.mkdir(cls.local_database_path)

        filename = os.path.join(cls.local_database_path, "DictyMutants")
        temp_file = os.path.join(cls.local_database_path, "DictyMutantsTemp")
        stream = urllib2.urlopen("http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID=all-mutants.txt")

        with open(temp_file, "wb") as file:
            shutil.copyfileobj(stream, file)

        if os.path.exists(filename):
            current = open(filename, "rb").read()
            temporary = open(temp_file, "rb").read()
            current.close()
            temporary.close()
            if current == temporary:
                os.remove(temp_file)
                return False

        os.rename(temp_file, filename)
        return True

    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            dicty = DictyMutants()
            cls._shared_dict = dicty.__dict__
        instance = DictyMutants.__new__(DictyMutants)
        instance.__dict__ = cls._shared_dict
        return instance

    def load(self, filename):
        file = open(filename, "rb")
        header = file.readline().rstrip()
        lines = file.read().splitlines()
        self._dicty_mutants = dict([(DictyMutant(line).id, DictyMutant(line)) for line in lines if line])

    def mutants(self):
        return self._dicty_mutants.values()

    def genes(self):
        return sorted(set(reduce(list.__add__, [self.mutant_genes(mutant.id) for mutant.id in self.mutants()], [])))

    def mutant_genes(self, mutant):
        return DictyMutant(self._dicty_mutants[mutant]).genes

    def gene_mutants(self):
        d = defaultdict(set)
        for mutant, genes in [(mutant, self.mutant_genes(mutant)) for mutant in self.mutants()]:
            for gene in genes:
                d[gene].add(mutant)
        return d

def mutants():
    """ Return all mutant objects
    """
    return DictyMutants.get_instance().mutants()

def genes():
    """ Return a set of all genes referenced in dictybase
    """
    return DictyMutants.get_instance().genes()

def mutant_genes(mutant):
    """ Return a set of all genes referenced by a mutant in dictybase
    """
    return DictyMutants.get_instance().mutant_genes(mutant)

def gene_mutants():
    """ Return a dictionary {gene: set(mutant_objects for mutant), ...}
    """
    return DictyMutants.get_instance().gene_mutants()

if  __name__  == "__main__":
    """
    Test whether the file contains only unique entries
    """
    entries = [ entry.id for entry in mutants() ]
    print len(set(entries)), len(entries)
    #print(mutants())
    #print(genes())
