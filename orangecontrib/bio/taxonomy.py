from __future__ import absolute_import, division

import os
import sys
import warnings

try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    from Orange.utils import environ
except ImportError:
    from orangecontrib.bio.utils import environ

from orangecontrib.bio.utils import serverfiles

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
    ("4577",   "Zea mays"),
    ("5476",   "Candida albicans")
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


def pickled_cache(filename=None, dependencies=[], version=1, maxSize=30,
                  pickleprotocol=pickle.HIGHEST_PROTOCOL):
    """
    Return a persistent cache function decorator.
    """
    def datetime_info(domain, filename):
        try:
            return serverfiles.info(domain, filename)["datetime"]
        except IOError:
            return serverfiles.ServerFiles().info(domain, filename)["datetime"]

    def cached(func):
        pytag = "py{0}.{1}".format(*sys.version_info[:2])
        if filename is None:
            cache_filename = os.path.join(
                environ.buffer_dir, func.__module__ + "_" + func.__name__ +
                "_" + pytag + "_cache.pickle")
        else:
            cache_filename = filename

        def f(*args, **kwargs):
            currentVersion = tuple([datetime_info(domain, file)
                                    for domain, file in dependencies] +
                                    [version, pytag])
            try:
                with open(cache_filename, "rb") as f:
                    cachedVersion, cache = pickle.load(f)
            except (OSError, IOError) as er:
                cachedVersion, cache = "no version", {}
            except ValueError as er:
                # unknown pickle format
                cachedVersion, cache = "no version", {}
            except Exception as er:
                warnings.warn(
                    "An error occurred while reading cache, using empty cache",
                    UserWarning)
                cachedVersion, cache = "no version", {}
            else:
                if cachedVersion != currentVersion:
                    cachedVersion, cache = "no version", {}

            allArgs = args + tuple([(key, tuple(value) if type(value) in [set, list] else value)\
                                     for key, value in kwargs.items()])
            if allArgs in cache:
                return cache[allArgs]
            else:
                res = func(*args, **kwargs)
                if len(cache) > maxSize:
                    del cache[next(iter(cache))]
                cache[allArgs] = res
                try:
                    with open(cache_filename, "wb") as f:
                        pickle.dump((currentVersion, cache), f,
                                    protocol=pickleprotocol)
                except (OSError, IOError) as err:
                    pass

                return res
        return f

    return cached


class Taxonomy(object):
    DOMAIN = "Taxonomy"
    FILENAME = "ncbi-taxonomy.sqlite"

    def __init__(self):
        from .ncbi.taxonomy import Taxonomy
        # Ensure the taxonomy db is downloaded.
        filename = serverfiles.localpath_download(self.DOMAIN, self.FILENAME)
        self._tax = Taxonomy(filename)

    def get_entry(self, id):
        try:
            return self._tax[id]
        except KeyError:
            raise UnknownSpeciesIdentifier(id)

    def search(self, string, onlySpecies=True, exact=False):
        res = self._tax.search(string, exact)
        if onlySpecies:
            res = [taxid for taxid in res
                   if self._tax[taxid].rank == "species"]
        return res

    def __iter__(self):
        return iter(self._tax)

    def __getitem__(self, id):
        return self.get_entry(id).name

    def other_names(self, id):
        return [(name, q) for name, q in self._tax[id].synonyms
                if q != "scientific name"]

    def rank(self, id):
        return self._tax[id].rank

    def parent(self, id):
        return self._tax[id].parent_tax_id

    def subnodes(self, id, levels=1):
        res = self._tax.child_tax_ids(id)
        if levels > 1:
            for child_id in res:
                res.extend(self.subnodes(child_id, levels - 1))
        return res

    def taxids(self):
        return list(self._tax)


@pickled_cache(None, [("Taxonomy", "ncbi-taxonomy.sqlite")], version=2)
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


@pickled_cache(None, [("Taxonomy", "ncbi-taxonomy.sqlite")], version=2)
def other_names(taxid):
    """
    Return a list of (name, name_type) tuples excluding the scientific name.

    Use :func:`name` to retrieve the scientific name.

    """
    return  Taxonomy().other_names(taxid)


@pickled_cache(None, [("Taxonomy", "ncbi-taxonomy.sqlite")], version=2)
def search(string, onlySpecies=True, exact=False):
    """ Search the NCBI taxonomy database for an organism.

    :param string: Search string.
    :param onlySpecies: Return only taxids of species (and subspecies).
    :param exact:  Return only taxids of organism that exactly match the string.
    """
    ids = Taxonomy().search(string, onlySpecies, exact)
    return list(ids)


@pickled_cache(None, [("Taxonomy", "ncbi-taxonomy.sqlite")], version=2)
def lineage(taxid):
    """ Return a list of taxids ordered from the topmost node (root) to taxid.
    """
    return Taxonomy()._tax.lineage(taxid)


def to_taxid(code, mapTo=None):
    """
    See if the code is a valid code in GO or KEGG and return a set of its taxids.
    """
    warnings.warn("'to_taxid' is deprecated", DeprecationWarning)

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
                elif r:
                    results.add(r)
            except Exception:
                pass

    if mapTo:
        mapTo = set(mapTo)
        if mapTo.intersection(results):
            return mapTo.intersection(results)

        mapped = set()
        for r in results:
            r_lin = lineage(r)
            if mapTo.intersection(r_lin):
                mapped.update(mapTo.intersection(r_lin))

        if not mapped:
            for m in mapTo:
                m_lin = lineage(m)
                if results.intersection(m_lin):
                    mapped.add(m)
        results = mapped

    return results


def taxids():
    """
    Returns a list of all (about half a million!) NCBI's taxonomy ID's.
    """
    return Taxonomy().taxids()


def ensure_downloaded(callback=None, verbose=True):
    """
    Retrieve the taxonomy database if not already downloaded.
    """
    serverfiles.localpath_download(
        Taxonomy.DOMAIN, Taxonomy.FILENAME,
        callback=callback, verbose=verbose)


if __name__ == "__main__":
    ids = search("Homo sapiens")
    print(ids)
    print(other_names(ids[0]))
