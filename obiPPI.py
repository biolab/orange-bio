"""obiPPI - Protein - protein interactions
"""

import os, sys
import xml.dom.minidom as minidom
import orngServerFiles

from obiKEGG import downloader

from collections import defaultdict

from obiTaxonomy import pickled_cache

        
class Interaction(object):
    def __init__(self, protein1, protein2, ref1=None, ref2=None, conf1=None, conf2=None):
        self.protein1, self.protein2 = protein1, protein2
        self.ref1, self.ref2 = ref1, ref2
        self.conf1, self.conf2 = conf1, conf2
        self.org1, self.org2 = None, None
    
class MIPS(object):
    VERSION = 1
    def __init__(self):
        self.load()
        
    def load(self):
        self.protein_names = defaultdict(set)
        self.refs = {}
        self.confidance = {}
        def process(element):
            d = {}
            participants = element.getElementsByTagName("proteinParticipant")
            proteins = []
            for protein in participants:
                interactor = protein.getElementsByTagName("proteinInteractor")[0]
                names = []
                for name in interactor.getElementsByTagName("shortLabel") + \
                            interactor.getElementsByTagName("fullName"):
                    names.append((name.tagName, name.childNodes[0].data))
                
                refs = []
                for ref in interactor.getElementsByTagName("primaryRef"):
                    refs += [(ref.tagName, ref.attributes.items())]
                org = dict(interactor.getElementsByTagName("organism")[0].attributes.items()).get("ncbiTaxId")
                conf = protein.getElementsByTagName("confidence")[0].attributes.items()
                proteins.append((names, refs, conf, org))
            interaction = Interaction(proteins[0][0][1][1], proteins[1][0][1][1])
            interaction.ref1, interaction.ref2 = proteins[0][1], proteins[1][1]
            interaction.conf1, interaction.conf2 = proteins[0][2], proteins[1][2]
            interaction.org1, interaction.org2 = proteins[0][3], proteins[1][3]
            
            self.protein_names[interaction.protein1].add(proteins[0][0][0][1])
            self.protein_names[interaction.protein2].add(proteins[1][0][0][1])
            
            return interaction 
            
        document = minidom.parse(orngServerFiles.localpath_download("PPI", "allppis.xml"))
        interactions = document.getElementsByTagName("interaction")
        self.interactions = [process(interaction) for interaction in interactions]
        
        self.protein_interactions = defaultdict(set)
        
        for inter in self.interactions:
            self.protein_names[inter.protein1] = dict(inter.ref1[0][1]).get("id")
            self.protein_names[inter.protein2] = dict(inter.ref2[0][1]).get("id")
            self.protein_interactions[inter.protein1].add(inter)
            self.protein_interactions[inter.protein2].add(inter) 
       
    def __iter__(self):
        return iter(self.interactions)
    
    @classmethod
    def download(cls):
        import urllib2, shutil
        src = urllib2.urlopen("http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz")
        dest = orngServerFiles.localpath("PPI", "mppi.gz")
        shutil.copyfileobj(src, open(dest, "wb"))
       
    @classmethod 
    @pickled_cache(None, [("PPI", "allppis.xml")], version=1)
    def _get_instance(cls):
        return MIPS()
    
    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_instance"):
            cls._instance= cls._get_instance()
        return cls._instance
    
def mips_interactions(protein = None):
    mips = MIPS.get_instance()
    if protein is None:
        return list(mips)
    else:
        return mips.protein_interactions.get(protein)

def mips_proteins():
    return set(MIPS.get_instance().protein_names.keys())

class BioGRIDInteraction(object):
    """ An object representing a BioGRID interaction. Each member of this object
    represents a data from a single column of BIOGRID-ALL.tab file.
    Attributes::
        - *interactor_a*    - BioGRID identifier
        - *interactor_b*    - BioGRID identifier
        - *official_symbol_a*    - An official symbol for *interactor_a*
        - *official_symbol_b*    - An official symbol for *interactor_b*
        - *aliases_for_a*    - Aliases separated by '|'
        - *aliases_for_b*    - Aliases separated by '|'
        - *experimental_system*     - Experimental system (see BioGRID documentation on www.thebiogrid.org for a list of valid entrys)
        - *source*    - 
        - *organism_a_id*    - NCBI Taxonomy identifier for *interactor_a*'s organism
        - *organism_b_id*    - NCBI Taxonomy identifier for *interactor_b*'s organism
    """
    __slots__ = ["interactor_a", "interactor_b", "official_symbol_a","official_symbol_b", "aliases_for_a", "aliases_for_b", "experimental_system", "source", "pubmed_id", "organism_a_id", "organism_b_id"]
    def __init__(self, line):
        for attr, val in zip(self.__slots__, line.split("\t")):
            setattr(self, attr, val)

class BioGRID(object):
    """ A BioGRID database interface
    Example::
        >>> ## finding all interactions for Homo sapiens sapiens
        >>> grid = BioGRID(case_insensitive=True)
        >>> proteins = proteins = biogrid.proteins() ## All proteins
        >>> proteins = [p for p in proteins if any(["9606" in [int.organism_a_id, int.organism_b_id] for int in grid.get(p)])]
    """
    VERSION = 1
    def __init__(self, case_insensitive=True):
        self.case_insensitive = case_insensitive
        self._case = (lambda name: name.lower()) if self.case_insensitive else (lambda name: name)
        self.load()
        
    def load(self):
        text = open(orngServerFiles.localpath_download("PPI", "BIOGRID-ALL.tab"), "rb").read()
        text = text.split("SOURCE\tPUBMED_ID\tORGANISM_A_ID\tORGANISM_B_ID\n", 1)[-1]
        self.interactions = [BioGRIDInteraction(line) for line in text.split("\n") if line.strip()]
        
        self.protein_interactions = defaultdict(set)
        self.protein_names = {}
        
        case = self._case

        def update(keys, value, collection):
            for k in keys:
                collection.setdefault(k, set()).add(value)
                
        for inter in self.interactions:
            update(map(case, [inter.official_symbol_a] + inter.aliases_for_a.split("|")), case(inter.interactor_a), self.protein_names)
            update(map(case, [inter.official_symbol_b] + inter.aliases_for_b.split("|")), case(inter.interactor_b), self.protein_names)
            
            self.protein_interactions[case(inter.interactor_a)].add(inter)
            self.protein_interactions[case(inter.interactor_b)].add(inter)
            
        self.protein_interactions = dict(self.protein_interactions)

        if case("N/A") in self.protein_names:
            del self.protein_names[case("N/A")]
        
    def proteins(self):
        """ Return all protein names in BioGRID (from INTERACTOR_A, and INTERACTOR_B columns) 
        """
        return self.protein_interactions.keys()
            
    def __iter__(self):
        """ Iterate over all BioGRIDInteraction objects
        """
        return iter(self.interactions)
    
    def __getitem__(self, key):
        """ Return a list of protein interactions that a protein is a part of 
        """
        key = self._case(key)
#        keys = self.protein_alias_matcher.match(key)
        if key not in self.protein_interactions:
            keys = self.protein_names.get(key, [])
        else:
            keys = [key]
        if keys:
            return list(reduce(set.union, [self.protein_interactions.get(k, []) for k in keys], set()))
        else:
            raise KeyError(key)
    
    def get(self, key, default=None):
        """ Return a list of protein interactions that a protein is a part of
        """
        key = self._case(key)
#        keys = self.protein_alias_matcher.match(key)
        if key not in self.protein_interactions:
            keys = self.protein_names.get(keys, [])
        else:
            keys = [key] 
        if keys:
            return list(reduce(set.union, [self.protein_interactions.get(k, []) for k in keys], set()))
        else:
            return default
        
    @classmethod
    def get_instance(cls):
        if getattr(cls, "_instance", None) is None:
            cls._instance = BioGRID()
        return cls._instance
    
def biogrid_interactions(name=None):
    """Return a list of protein interactions (BioGRIDInteraction objects) that a protein is a part of
    """ 
    if name:
        return list(BioGRID.get_instance().get(name, set()))
    else:
        return BioGRID.get_instance().interactions
    
def biogrid_proteins():
    """ Return all protein names in BioGRID (from INTERACTOR_A, and INTERACTOR_B columns)
    """
    return BioGRID.get_instance().proteins()


if __name__ == "__main__":
    for protein in mips_proteins():
        print "Protein", protein, "interacts with", 
        print ",".join(set(reduce(list.__add__, [[inter.protein1, inter.protein2] for inter in mips_interactions(protein)], [])) -set([protein]))
            