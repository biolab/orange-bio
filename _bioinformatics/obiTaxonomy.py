from __future__ import absolute_import, division

from collections import defaultdict
import cPickle, os, shutil, sys, StringIO, tarfile, urllib2

from Orange.orng import orngEnviron, orngServerFiles

from . import obiData, obiGenomicsUpdate

# list of common organisms from http://www.ncbi.nlm.nih.gov/Taxonomy
def common_taxids():
    """Return taxonomy IDs for common organisms."""
    return ["3702",  # Arabidopsis thaliana
            "9913",  # Bos taurus
            "6239",  # Caenorhabditis elegans
            "3055",  # Chlamydomonas reinhardtii
            "7955",  # Danio rerio (zebrafish)
            "352472", # Dictyostelium discoideum
            "7227",  # Drosophila melanogaster
            "562",   # Escherichia coli
            "11103", # Hepatitis C virus
            "9606",  # Homo sapiens
            "10090", # Mus musculus
            "2104",  # Mycoplasma pneumoniae
            "4530",  # Oryza sativa
            "5833",  # Plasmodium falciparum
            "4754",  # Pneumocystis carinii
            "10116", # Rattus norvegicus
            "4932",  # Saccharomyces cerevisiae
            "4896",  # Schizosaccharomyces pombe
            "31033", # Takifugu rubripes
            "8355",  # Xenopus laevis
            "4577",  # Zea mays
             ] 

def essential_taxids():
    """Return taxonomy IDs for organisms that are included in (default) Orange Bioinformatics installation."""
    return ["352472", # Dictyostelium discoideum
            "7227",  # Drosophila melanogaster
            "9606",  # Homo sapiens
            "10090", # Mus musculus
            "4932",  # Saccharomyces cerevisiae
            ] 

default_database_path = os.path.join(orngEnviron.bufferDir, "bigfiles", "Taxonomy")

class MultipleSpeciesException(Exception):
    pass

class UnknownSpeciesIdentifier(Exception):
    pass

def pickled_cache(filename=None, dependencies=[], version=1, maxSize=30):
    """ Return a cache function decorator. 
    """
    def datetime_info(domain, filename):
        try:
            return orngServerFiles.info(domain, filename)["datetime"]
        except IOError:
            return orngServerFiles.ServerFiles().info(domain, filename)["datetime"]
            
    def cached(func):
        default_filename = os.path.join(orngEnviron.bufferDir, func.__module__ + "_" + func.__name__ + "_cache.pickle")
        def f(*args, **kwargs):
            currentVersion = tuple([datetime_info(domain, file) for domain, file in dependencies]) + (version,)
            try:
                cachedVersion, cache = cPickle.load(open(filename or default_filename, "rb"))
                if cachedVersion != currentVersion:
                    cache = {}
            except IOError, er:
                cacheVersion, cache = "no version", {}
            allArgs = args + tuple([(key, tuple(value) if type(value) in [set, list] else value)\
                                     for key, value in kwargs.items()])
            if allArgs in cache:
                return cache[allArgs]
            else:
                res = func(*args, **kwargs)
                if len(cache) > maxSize:
                    del cache[iter(cache).next()]
                cache[allArgs] = res
                cPickle.dump((currentVersion, cache), open(filename or default_filename, "wb"), protocol=cPickle.HIGHEST_PROTOCOL)
                return res
        return f

    return cached

def cached(func):
    """Cached one arg method
    """
    def f(self, arg):
        if arg not in self._cache:
            self._cache[arg] = func(self, arg)
        return self._cache[arg]
    f._cache = {}
    f.__name__ = "Cached " + func.__name__
    return f
    
class TextDB(object):
    entry_start_string = chr(255)
    entry_end_string = chr(254)+"\n"
    entry_separator_string = chr(253)
    
    @property
    def _text_lower(self):
        if id(self._text) == self._lower_text_id:
            return self._lower_text
        else:
            self._lower_text_id = id(self._text)
            self._lower_text = self._text.lower()
            return self._lower_text
        
    def __init__(self, file=None, **kwargs):
        self._text = ""
        self._lower_text_id = id(self._text) - 1
        self._cache = {}
        self.__dict__.update(kwargs)
        
        if file != None:
            self._text = open(file, "rb").read()

    def _find_all(self, string, start=0, text=None, unique=True):
        text = text if text != None else self._text_lower
        while True:
            index = text.find(string, start)
            if index != -1:
                yield index
                if unique:
                    start = text.find(self.entry_start_string, index + 1)
                else:
                    start = index + 1
            else:
                raise StopIteration

    def _get_entry_at(self, index, text=None):
        text = text if text != None else self._text
        start = text.rfind(self.entry_start_string, 0, index + 1)
        end = text.find(self.entry_end_string, index)
        return self._text[start+1:end]

    @cached
    def get_entry(self, id):
        try:
            index = self._find_all(self.entry_start_string + id + self.entry_separator_string).next()
        except StopIteration:
            raise KeyError, id
        return self._get_entry_at(index)
                
    def search(self, string):
        string = string.lower()
        res = []
        for idx in self._find_all(string):
            entry = self._get_entry_at(idx)
            id , rest = entry.split(self.entry_separator_string, 1)
            self._cache[id] = entry
            res.append(id)
        return res

    def insert(self, entry):
        self._text += self.entry_start_string + self.entry_separator_string.join(entry) + self.entry_end_string

    def __iter__(self):
        for idx in self._find_all(self.entry_start_string):
            entry = self._get_entry_at(idx)
            if entry:
                yield entry.split(self.entry_separator_string ,1)[0]

    def __getitem__(self, id):
        entry = self.get_entry(id)
        return entry.split(self.entry_separator_string)[1:]

    def __setitem__(self, id, entry):
        self.insert([id] + list(entry))

    def create(self, filename):
        f = open(filename, "wb")
        f.write(self._text)
        def write(entry):
            f.write(self.entry_start_string + self.entry_separator_string.join(entry) + self.entry_end_string)
        return write
    
class Taxonomy(object):
    __shared_state = {"_text":None, "_info":None}
    def __init__(self):
        self.__dict__ = self.__shared_state
        if not self._text:
            self.Load()
            
    def Load(self):
        try:
            self._text = TextDB(os.path.join(default_database_path, "ncbi_taxonomy.tar.gz", "ncbi_taxonomy.db"))
            self._info = TextDB(os.path.join(default_database_path, "ncbi_taxonomy.tar.gz", "ncbi_taxonomy_inf.db"))
            return
        except Exception, ex:
            pass
        try:
            from . import orngServerFiles as sf
            sf.download("Taxonomy", "ncbi_taxonomy.tar.gz")
            self._text = TextDB(os.path.join(default_database_path, "ncbi_taxonomy.tar.gz", "ncbi_taxonomy.db"))
            self._info = TextDB(os.path.join(default_database_path, "ncbi_taxonomy.tar.gz", "ncbi_taxonomy_inf.db"))
            return
        except Exception, ex:
            raise

    def get_entry(self, id):
        try:
            entry = self._text[id]
        except KeyError:
            raise UnknownSpeciesIdentifier
        return entry
                
    def search(self, string, onlySpecies=True):
        res = self._text.search(string)
        if onlySpecies:
            res = [r for r in res if "species" in self._text[r][1]]
        return res

    def __iter__(self):
        return iter(self._text)

    def __getitem__(self, id):
        entry = self.get_entry(id)
        return entry[2] ## item with index 2 is allways scientific name

    def other_names(self, id):
        entry = self.get_entry(id)
        info = self._info[id]
        names = entry[2:] ## index 2 and larger are names
        return list(zip(names, info))[1:] ## exclude scientific name

    def rank(self, id):
        entry = self.get_entry(id)
        return entry[1]

    def parent(self, id):
        entry = self.get_entry(id)
        return entry[0]

    def subnodes(self, id, levels=1):
        res = self._text.search(self._text.entry_separator_string + id + self._text.entry_separator_string)
        res = [r for r in res if self.get_entry(r)[0] == id]
        if levels > 1:
            for r in list(res):
                res.extend(self.subnodes(r, levels-1))
        return res

    def taxids(self):
        return list(self)
    
    @staticmethod
    def ParseTaxdumpFile(file=None, outputdir=None, callback=None):
        from cStringIO import StringIO
        if file == None:
            file = tarfile.open(None, "r:gz", StringIO(urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz").read()))
        if type(file) == str:
            file = tarfile.open(file)
        names = file.extractfile("names.dmp").readlines()
        nodes = file.extractfile("nodes.dmp").readlines()
        namesDict = defaultdict(list)
        for line in names:
            if not line.strip():
                continue
            line = line.rstrip("\t\n|").split("\t|\t")
            id, name, unique_name, name_class = line
            if unique_name:
                namesDict[id].append((unique_name , name_class))
            else:
                namesDict[id].append((name , name_class))

        nodesDict = {}
        for line in nodes:
            if not line.strip():
                continue
            line = line.split("\t|\t")[:3]
            id, parent, rank = line
            nodesDict[id] = (parent, rank)
        
        name_class_codes = defaultdict(iter(range(255)).next)
        name_class_codes["scientific name"] ## Force scientific name to be first
        if outputdir == None:
            outputdir = default_database_path
        text = TextDB().create(os.path.join(outputdir, "ncbi_taxonomy.db"))
        info = TextDB().create(os.path.join(outputdir, "ncbi_taxonomy_inf.db"))
        milestones = set(range(0, len(namesDict), max(len(namesDict)/100, 1)))
        for i, (id, names) in enumerate(namesDict.items()):
            parent, rank = nodesDict[id]
            ## id, parent and rank go first
            entry = [id, parent, rank]
            ## all names and name class codes pairs follow ordered so scientific name is first
            names = sorted(names, key=lambda (name, class_): name_class_codes[class_])
            entry.extend([name for name ,class_ in names])
            info_entry = [id] + [class_ for name, class_ in names]
            text(entry)
            info(info_entry)
            if callback and i in milestones:
                callback(i)
    
def name(taxid):
    """ Return the scientific name for organism with taxid.
    """
    return Taxonomy()[taxid]

def other_names(taxid):
    """ Return a list of (name, name_type) tuples but exclude the scientific name.
    """
    return  Taxonomy().other_names(taxid)

@pickled_cache(None, [("Taxonomy", "ncbi_taxonomy.tar.gz")], version=1)
def search(string, onlySpecies=True, exact=False):
    """ Search the NCBI taxonomy database for an organism
    Arguments::
            - *string*      Search string
            - *onlySpecies* Return only taxids of species (and subspecies)
            - *exact*       Return only taxids of organism that exactly match the string
    """
    ids = Taxonomy().search(string, onlySpecies)
    if exact:
        ids = [id for id in ids if string in [name(id)] + [t[0] for t in other_names(id)]]
    return ids

def lineage(taxid):
    """ Return a list of taxids ordered from the topmost node (root) to taxid.
    """
    tax = Taxonomy()
    result = [taxid]
    while True:
        parent = tax.parent(result[-1])
        result.append(parent)
        if tax[parent] == "root" or parent=="1":
            break
    result.reverse()
    return result
    
def to_taxid(code, mapTo=None):
    """ See if the code is a valid code in any database and return a set of its taxids.
    """
    from . import obiKEGG, obiGO
    results = set()
    for test in [obiKEGG.to_taxid, obiGO.to_taxid]:
        try:
            r = test(code)
            if type(r) == set:
                results.update(r)
            else:
                results.add(r)
        except Exception, ex:
            pass

    if mapTo:
        mapped = [[parent_id for parent_id in mapTo if parent_id in lineage(r)].pop() for r in results if any(parent_id in lineage(r) for parent_id in mapTo)]
        if not mapped:
            t = Taxonomy()
            subnodes = dict([(r, t.subnodes(r, 4)) for r in results])
            mapped = [id for id in mapTo if any(id in subnodes[r] for r in results)]
        results = mapped

    return results

def taxids():
    return Taxonomy().taxids()

from . import obiGenomicsUpdate

class Update(obiGenomicsUpdate.Update):
    def GetDownloadable(self):
        return [Update.UpdateTaxonomy]
    
    def IsUpdatable(self, func, args):
        from datetime import datetime
        from . import obiData
        if func == Update.UpdateTaxonomy:
##            stream = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
##            date = datetime.strptime(stream.headers.get("Last-Modified"), "%a, %d %b %Y %H:%M:%S %Z")
            ftp = obiData.FtpWorker("ftp.ncbi.nih.gov")
            size, date = ftp.statFtp("pub/taxonomy/taxdump.tar.gz")
            return date > self.GetLastUpdateTime(func, args)

    def UpdateTaxonomy(self):
        Taxonomy.ParseTaxdumpFile(outputdir=self.local_database_path)
        import tarfile
        tFile = tarfile.open(os.path.join(self.local_database_path, "ncbi_taxonomy.tar.gz"), "w:gz")
        tFile.add(os.path.join(self.local_database_path, "ncbi_taxonomy.db"), "ncbi_taxonomy.db")
        tFile.add(os.path.join(self.local_database_path, "ncbi_taxonomy_inf.db"), "ncbi_taxonomy_inf.db")
        tFile.close()

if __name__ == "__main__":
    ids = search("Homo sapiens")
    print ids
    print other_names(ids[0])
    
