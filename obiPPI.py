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


if __name__ == "__main__":
    for protein in mips_proteins():
        print "Protein", protein, "interacts with", 
        print ",".join(set(reduce(list.__add__, [[inter.protein1, inter.protein2] for inter in mips_interactions(protein)], [])) -set([protein]))
            