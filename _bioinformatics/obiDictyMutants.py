import os
import urllib2
import shutil

from collections import defaultdict

from Orange.orng import orngServerFiles

class DictyMutant(object):

    def __init__(self, mutant_line):
        mutant = mutant_line.split("\t")
        self.name = mutant[0]
        self.descriptor = set(mutant[1].split("/"))
        self.genes = mutant[2].split(" | ")
        self.phenotypes = mutant[3].split(" | ")
        self.null = False
        self.overexp = False
        self.multiple = False
        self.develop = False
        self.other = False
 
class DictyMutants(object):
    VERSION=1
    DEFAULT_DATABASE_PATH = orngServerFiles.localpath("DictyMutants") #use a default folder for storing the genesets

    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path is not None else self.DEFAULT_DATABASE_PATH

        if not os.path.exists(self.local_database_path):
            os.mkdir(self.local_database_path)
            
        _mutants = self.prepare_mutants()
 
    def update_file(self, name):
        url = "http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID="
        filename = os.path.join(self.local_database_path, name)
        temp_file = os.path.join(self.local_database_path, name + "_temp")
        stream = urllib2.urlopen(url + name)
    
        with open(temp_file, "wb") as file:
            shutil.copyfileobj(stream, file)
    
        os.rename(temp_file, filename)
        return filename
    
    def load_mutants(self, file):
        data = open(file)
        data_header = data.readline()
        data = data.read()
        return data.splitlines()
       
    def prepare_mutants(self):   
        all_mutants = self.load_mutants(self.update_file("all-mutants.txt"))
        null_mutants = self.load_mutants(self.update_file("null-mutants.txt"))
        overexp_mutants = self.load_mutants(self.update_file("overexpression-mutants.txt"))
        multiple_mutants = self.load_mutants(self.update_file("multiple-mutants.txt"))
        develop_mutants = self.load_mutants(self.update_file("developmental-mutants.txt"))
        other_mutants = self.load_mutants(self.update_file("other-mutants.txt"))
   
        _mutants = [DictyMutant(mutant) for mutant in all_mutants]
        
        the_nulls = set([DictyMutant(line).name for line in null_mutants])
        the_overexps = set([DictyMutant(line).name for line in overexp_mutants])
        the_multiples = set([DictyMutant(line).name for line in multiple_mutants])
        the_develops = set([DictyMutant(line).name for line in develop_mutants])
        the_others = set([DictyMutant(line).name for line in other_mutants])

        for mutant in _mutants:
            if mutant.name in the_nulls: mutant.null = True
            if mutant.name in the_overexps: mutant.overexp = True 
            if mutant.name in the_multiples: mutant.multiple = True
            if mutant.name in the_develops: mutant.develop = True
            if mutant.name in the_others: mutant.other = True
       
        self._mutants = {x: x for x in _mutants}

    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            dicty = DictyMutants()
            cls._shared_dict = dicty.__dict__
        instance = DictyMutants.__new__(DictyMutants)
        instance.__dict__ = cls._shared_dict
        return instance

    def mutants(self):
        return self._mutants.keys()

    def genes(self):
        return sorted(set(reduce(list.__add__, [self.mutant_genes(mutant) for mutant in self.mutants()], [])))

    def mutant_genes(self, mutant):
        return self._mutants[mutant].genes

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
    print(mutants())    
