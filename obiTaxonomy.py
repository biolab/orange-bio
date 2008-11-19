import urllib2
import os, sys, shutil
import cPickle
import tarfile
import StringIO
import obiGenomicsUpdate
import obiData
import orngEnviron

from collections import defaultdict

default_database_path = os.path.join(orngEnviron.bufferDir, "bigfiles", "Taxonomy")

class MultipleSpeciesException(Exception):
    pass

class UnknownSpeciesIdentifier(Exception):
    pass

class Taxonomy(object):
    __instance = None
    def __init__(self, file=None):
        if file:
            self.ParseFile(file)

    def ParseFile(self, file):
        """Parse the taxdump.tar.gz file from NCBI
        """
        if type(file) == str:
            tfile = tarfile.open(file)
        names = tfile.extractfile("names.dmp").readlines()
        nodes = tfile.extractfile("nodes.dmp").readlines()
        self.names = namesDict = defaultdict(list)
        for line in names:
            if not line.strip():
                continue
            line = line.rstrip("\t\n|").split("\t|\t")
            id, name, unique_name, name_class = line
            
            if unique_name:
                namesDict[id].append((unique_name , name_class))
            else:
                namesDict[id].append((name , name_class))

        self.nodes = nodesDict = {}
        for line in nodes:
            if not line.strip():
                continue
            line = line.split("\t|\t")[:3]
            id, parent, rank = line
            nodesDict[id] = (parent, rank)

    @classmethod
    def Load(cls):
        if not Taxonomy.__instance:
            Taxonomy.__instance = cPickle.load(open(os.path.join(default_database_path, "taxonomy.pickle")))
        return Taxonomy.__instance

    @staticmethod
    def DownloadTaxonomy(file=None, progressCallback=None):
        if file == None:
            file = os.path.join(default_database_path, "taxdump.tar.gz")
        if type(file) == str:
            file = open(file, "wb")
        fd = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
        shutil.copyfileobj(fd, file)
        file.close()
        tax = Taxonomy(os.path.join(default_database_path, "taxdump.tar.gz"))
        cPickle.dump(tax, open(os.path.join(default_database_path, "taxonomy.pickle"), "wb"))
    
    def __iter__(self):
        return self.names.__iter__()

    def __getitem__(self, id):
        return self.names[id]
    
def name(taxid):
    tax = Taxonomy.Load()
    names = [name for name, type in tax.names[taxid] if type == "scientific name"]
    return names[0]

def other_names(taxid):
    tax = Taxonomy.Load()
    return tax.names[taxid]

def search(string, onlySpecies=True):
    tax = Taxonomy.Load()
    string = string.lower()
    match = []
    for id, names in tax.names.iteritems():
        if any(string in name.lower() for name, type in names):
            if onlySpecies == False or tax.nodes[id][1] == "species":
                match.append(id)
    return match

if __name__ == "__main__":
    ids = search("Homo sapiens")
    print ids
    print other_names(ids[0])
    