import urllib2
import os, sys, shutil
import tarfile
import StringIO
import obiGenomicsUpdate
import obiData
import orngEnviron

from collections import defaultdict

default_database_path = os.path.join(orngEnviron.bufferDir, "bigfiles", "Taxonomy")

class MultipleSpeciesExceptions(Exception):
    pass

def _download_taxonomy_NCBI(file=None, progressCallback=None):
    if file == None:
        file = os.path.join(default_database_path, "taxdump.tar.gz")
    if type(file) == str:
        file = tarfile.open(file, "w:gz")
    fd = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
    buffer = StringIO.StringIO()
    shutil.copyfileobj(fd, buffer)
    buffer.seek(0)
    tFile = tarfile.open(None, "r:gz", buffer)
    namesInfo = tFile.getmember("names.dmp")
    nodesInfo = tFile.getmember("nodes.dmp")
    file.addfile(namesInfo, tFile.extractfile(namesInfo))
    file.addfile(nodesInfo, tFile.extractfile(nodesInfo))
    file.close()

class _taxonomy(object):
    def __init__(self, names, nodes):
        self.names = names
        self.nodes = nodes
        
def get_taxonomy(_cached=_taxonomy(None, None)):
    if _cached.names and _cached.nodes:
        return _cached
    filename = os.path.join(default_database_path, "taxdump.tar.gz")
    if not os.path.isfile(filename):
        try:
            os.mkdir(os.path.dirname(filename))
        except Exception:
            pass
        _download_taxonomy_NCBI(tarfile.open(os.path.join(default_database_path, "taxdump.tar.gz"), "w"))
    tfile = tarfile.open(filename)
    names = tfile.extractfile("names.dmp").read()
    nodes = tfile.extractfile("nodes.dmp").read()
    namesDict = defaultdict(list)
    for line in names:
        if not line:
            continue
        line = line.rstrip("\t\n|").split("\t|\t")
        id, name, unique_name, name_class = line
        if unique_name:
            namesDict[id].append((unique_name , name_class))
        else:
            namesDict[id].append((name , name_class))

    nodesDict = {}
    for line in nodes:
        if not line:
            continue
        line = line.split("\t|\t")[:3]
        id, parent, rank = line
        nodesDict[id] = (parent, rank)
        
    _cached.names = dict(namesDict)
    _cached.nodes = nodesDict
    return _cached
    
def name(taxid):
    tax = get_taxonomy()
    names = [name for name, type in tax.names[taxid] if type == "scientific name"]
    return names[0]

def other_names(taxid):
    tax = get_taxonomy()
    return tax.names[taxid]

def search(string, onlySpecies=True):
    tax = get_taxonomy()
    string = string.lower()
    match = []
    for id, names in tax.names.items():
        if any(string in name.lower() for name, type in names):
            if onlySpecies == False or tax.nodes[id][1] == "species":
                match.append(id)
    return match

if __name__ == "__main__":
    ids = search("Homo sapiens")
    print ids
    print other_names(ids[0])
    