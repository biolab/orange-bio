import os
import urllib2
import shutil

from collections import defaultdict

from Orange.orng import orngServerFiles

class dicty_mutant(object):

    def __init__(self, mutant_line):
        mutant = mutant_line.split("\t")
        self.id = mutant[0]
        self.descriptor = set(mutant[1].split("/"))
        self.genes = mutant[2].split(" | ")
        self.phenotypes = mutant[3].split(" | ")


class dicty_mutants(object):
    VERSION=1
    DEFAULT_DATABASE_PATH = orngServerFiles.localpath("DictyMutants") #use a default folder for storing the genesets

    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path is not None else self.DEFAULT_DATABASE_PATH

        if not os.path.exists(self.local_database_path):
            download_from_dictybase(self, self.local_database_path)

        filename = os.path.join(cls.local_database_path, "DictyMutants")
        self.load(filename)

    @classmethod
    def download_from_dictybase(cls, local_database_path=None):
        cls.local_database_path = local_database_path if local_database_path is not None else cls.DEFAULT_DATABASE_PATH

        if not os.path.exists(cls.local_database_path):
            os.mkdir(cls.local_database_path)

        filename = os.path.join(cls.local_database_path, "DictyMutants")
        temp_file = os.path.join(cls.local_database_path, "DictyMutantsTemp")
        stream = urllib2.urlopen("http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID=all-mutants.txt")

        if not os.path.exists(filename):
            with open(filename, "wb") as file:
                shutil.copyfileobj(stream, file)
        else:
            toremove = False
            with open(temp_file, "wb") as file:
                shutil.copyfileobj(stream, file)
            current = open(filename, "rb")
            temporary = open(temp_file, "rb")
            if current.read() == temporary.read():
                toremove = True
            current.close()
            temporary.close()

            if toremove:
                os.remove(temp_file)
            else:
                os.rename(temp_file, filename)

    def load(self, filename):
        file = open(filename, "rb")
        header = file.readline().rstrip()
        lines = file.read().splitlines()
        self._dicty_mutants = dict([(dicty_mutant(line), line) for line in lines if line])

    def mutants(self):
        return self._dicty_mutants.keys()

    def genes(self):
        return sorted(set(reduce(list.__add__, [self.mutant_genes(mutant) for mutant in self.mutants()], [])))

    def mutant_genes(self, mutant):
        return dicty_mutant(self._dicty_mutants[mutant]).genes

    def gene_mutants(self):
        d = defaultdict(set)
        for mutant, genes in [(mutant, self.mutant_genes(mutant)) for mutant in self.mutants()]:
            for gene in genes:
                d[gene].add(mutant)
        return d

def mutants():
    """ Return all mutant objects
    """
    return dicty_mutants.get_instance().mutants()

def genes():
    """ Return a set of all genes referenced in dictybase
    """
    return dicty_mutants.get_instance().genes()

def mutant_genes(mutant):
    """ Return a set of all genes referenced by a mutant in dictybase
    """
    return OMIM.get_instance().mutant_genes(mutant)

def gene_mutants():
    """ Return a dictionary {gene: set(disease_objects for gene), ...}
    """
    return OMIM.get_instance().gene_diseases()


