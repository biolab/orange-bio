from __future__ import absolute_import, division

from collections import defaultdict

try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import  os, shutil, sys, tarfile

try:
    from Orange.utils import environ
except ImportError:
    from .utils import environ

from .utils import serverfiles

import six

_COMMON_NAMES = (
    ("3702",   "Arabidopsis thaliana"),
    ("9913",   "Bos taurus"),
    ("6239",   "Caenorhabditis elegans"),
    ("3055",   "Chlamydomonas reinhardtii"),
    ("7955",   "Danio rerio"),
    ("352472", "Dictyostelium discoideum AX4"),
    ("7227",   "Drosophila melanogaster"),
    ("562",    "Escherichia coli"),
    ("11103",  "Hepatitis C virus"),
    ("9606",   "Homo sapiens"),
    ("10090",  "Mus musculus"),
    ("2104",   "Mycoplasma pneumoniae"),
    ("4530",   "Oryza sativa"),
    ("5833",   "Plasmodium falciparum"),
    ("4754",   "Pneumocystis carinii"),
    ("10116",  "Rattus norvegicus"),
    ("4932",   "Saccharomyces cerevisiae"),
    ("4896",   "Schizosaccharomyces pombe"),
    ("31033",  "Takifugu rubripes"),
    ("8355",   "Xenopus laevis"),
    ("4577",   "Zea mays")
)

_COMMON_NAMES_MAPPING = dict(_COMMON_NAMES)


# list of common organisms from http://www.ncbi.nlm.nih.gov/Taxonomy
def common_taxids():
    """Return taxonomy IDs for common organisms."""
    # Sorted lexicographically by names
    return [taxid for taxid, _ in _COMMON_NAMES]


def common_taxid_to_name(taxid):
    """Return a name for a common organism taxonomy id."""
    return _COMMON_NAMES_MAPPING[taxid]


def taxname_to_taxid(name):
    """Return taxonomy ID for a taxonomy name."""
    name_to_taxid = dict(map(reversed, _COMMON_NAMES))

    if name in name_to_taxid:
        return name_to_taxid[name]
    return None


def essential_taxids():
    """Return taxonomy IDs for organisms that are included in (default)
    Orange Bioinformatics installation.
    """
    return ["352472", # Dictyostelium discoideum
            "7227",  # Drosophila melanogaster
            "9606",  # Homo sapiens
            "10090", # Mus musculus
            "4932",  # Saccharomyces cerevisiae
            ]


def shortname(taxid):
    """ Short names for common_taxids organisms """
    names = {
    "3702"  : ["arabidopsis", "thaliana", "plant"],
    "9913"  : ["cattle", "cow"],
    "6239"  : ["nematode", "roundworm"],
    "3055"  : ["algae"],
    "7955"  : ["zebrafish"],
    "352472": ["dicty", "amoeba", "slime mold"],
    "7227"  : ["fly", "fruit fly", "vinegar fly"],
    "562"   : ["ecoli", "coli", "bacterium"],
    "11103" : ["virus, hepatitis"],
    "9606"  : ["human"],
    "10090" : ["mouse", "mus"],
    "2104"  : ["bacterium", "mycoplasma"],
    "4530"  : ["asian rice", "rice", "cereal", "plant"],
    "5833"  : ["plasmodium", "malaria", "parasite"],
    "4754"  : ["pneumonia", "fungus"],
    "10116" : ["rat", "laboratory rat"],
    "4932"  : ["yeast", "baker yeast", "brewer yeast"],
    "4896"  : ["yeast", "fission yeast"],
    "31033" : ["fish", "pufferfish"],
    "8355"  : ["frog", "african clawed frog"],
    "4577"  : ["corn", "cereal grain", "plant"]
    }
    if taxid in names:
        return names[taxid]
    return []

class MultipleSpeciesException(Exception):
    pass

class UnknownSpeciesIdentifier(Exception):
    pass

def pickled_cache(filename=None, dependencies=[], version=1, maxSize=30):
    """ Return a cache function decorator. 
    """
    def datetime_info(domain, filename):
        try:
            return serverfiles.info(domain, filename)["datetime"]
        except IOError:
            return serverfiles.ServerFiles().info(domain, filename)["datetime"]
            
    def cached(func):
        default_filename = os.path.join(environ.buffer_dir, func.__module__ + "_" + func.__name__ + "_cache.pickle")
        def f(*args, **kwargs):
            currentVersion = tuple([datetime_info(domain, file) for domain, file in dependencies]) + (version,)
            try:
                cachedVersion, cache = pickle.load(open(filename or default_filename, "rb"))
                if cachedVersion != currentVersion:
                    cache = {}
            except IOError as er:
                cacheVersion, cache = "no version", {}
            allArgs = args + tuple([(key, tuple(value) if type(value) in [set, list] else value)\
                                     for key, value in kwargs.items()])
            if allArgs in cache:
                return cache[allArgs]
            else:
                res = func(*args, **kwargs)
                if len(cache) > maxSize:
                    del cache[next(iter(cache))]
                cache[allArgs] = res
                pickle.dump((currentVersion, cache), open(filename or default_filename, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
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
            if six.PY3:
                self._text = open(file, "rt", encoding="latin1").read()
            else:
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
            index = next(self._find_all(self.entry_start_string + id + self.entry_separator_string))
        except StopIteration:
            raise KeyError(id)
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
        path = serverfiles.localpath_download("Taxonomy", "ncbi_taxonomy.tar.gz")
        self._text = TextDB(os.path.join(path, "ncbi_taxonomy.db"))
        self._info = TextDB(os.path.join(path, "ncbi_taxonomy_inf.db"))

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
        import Orange.utils
        if file == None:
            so = StringIO()
            Orange.utils.wget("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", dst_obj=so)
            file = tarfile.open(None, "r:gz", StringIO(so.getvalue()))
            so.close()
        elif type(file) == str:
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
        
        if outputdir == None:
            outputdir = serverfiles.localpath("Taxonomy")
        text = TextDB().create(os.path.join(outputdir, "ncbi_taxonomy.db"))
        info = TextDB().create(os.path.join(outputdir, "ncbi_taxonomy_inf.db"))
        milestones = set(range(0, len(namesDict), max(int(len(namesDict)/100), 1)))
        for i, (id, names) in enumerate(namesDict.items()):
            parent, rank = nodesDict[id]
            ## id, parent and rank go first
            entry = [id, parent, rank]
            ## all names and name class codes pairs follow ordered so scientific name is first
            names = sorted(names, key=lambda x: (not x[1] == "scientific name", x[1], x[0]))
            print(names)
            entry.extend([name for name ,class_ in names])
            info_entry = [id] + [class_ for name, class_ in names]
            text(entry)
            info(info_entry)
            if callback and i in milestones:
                callback(i)


@pickled_cache(None, [("Taxonomy", "ncbi_taxonomy.tar.gz")], version=1)
def name(taxid):
    """
    Return the scientific name for organism with taxid.
    """
    # Most of the lookups will be for the common names, so in most
    # situations we can avoid loading the taxonomy.
    if taxid in _COMMON_NAMES:
        return _COMMON_NAMES[taxid]
    else:
        return Taxonomy()[taxid]


@pickled_cache(None, [("Taxonomy", "ncbi_taxonomy.tar.gz")], version=1)
def other_names(taxid):
    """
    Return a list of (name, name_type) tuples excluding the scientific name.

    Use :func:`name` to retrieve the scientific name.

    """
    return  Taxonomy().other_names(taxid)


@pickled_cache(None, [("Taxonomy", "ncbi_taxonomy.tar.gz")], version=1)
def search(string, onlySpecies=True, exact=False):
    """ Search the NCBI taxonomy database for an organism.

    :param string: Search string.
    :param onlySpecies: Return only taxids of species (and subspecies).
    :param exact:  Return only taxids of organism that exactly match the string.
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
    """
    See if the code is a valid code in GO or KEGG and return a set of its taxids.
    """
    try:
        name(code)
        results = set([code])
    except UnknownSpeciesIdentifier:
        results = set()
        from . import kegg, go
        for test in [kegg.to_taxid, go.to_taxid]:
            try:
                r = test(code)
                if type(r) == set:
                    results.update(r)
                else:
                    results.add(r)
            except Exception as ex:
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
    """ Returns a list of all (about half a million!) NCBI's taxonomy ID's. """
    return Taxonomy().taxids()


def ensure_downloaded(callback=None, verbose=True):
    """
    Retrieve the taxonomy database if not already downloaded.
    """
    serverfiles.localpath_download("Taxonomy", "ncbi_taxonomy.tar.gz",
                                   callback=callback, verbose=verbose)


if __name__ == "__main__":
    ids = search("Homo sapiens")
    print(ids)
    print(other_names(ids[0]))
