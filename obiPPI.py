"""obiPPI - Protein - protein interactions
"""

import os, sys
import xml.dom.minidom as minidom
import orngServerFiles

from obiKEGG import downloader

from collections import defaultdict

        
class Interaction(object):
    def __init__(self, protein1, protein2, ref1=None, ref2=None, conf1=None, conf2=None):
        self.protein1, self.protein2 = protein1, protein2
        self.ref1, self.ref2 = ref1, ref2
        self.conf1, self.conf2 = conf1, conf2
    
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
                conf = protein.getElementsByTagName("confidence")[0].attributes.items()
                proteins.append((names, refs, conf))
            interaction = Interaction(proteins[0][0][1][1], proteins[1][0][1][1])
            interaction.ref1, interaction.ref2 = proteins[0][1], proteins[1][1]
            interaction.conf1, interaction.conf2 = proteins[0][2], proteins[1][2]
            
            self.protein_names[interaction.protein1].add(proteins[0][0][0][1])
            self.protein_names[interaction.protein2].add(proteins[1][0][0][1])
            
            self.protein_names[interaction.protein1].add(proteins[0][0][0][1])
            self.protein_names[interaction.protein1].add(proteins[0][0][0][1])
            
            return interaction 
            
        document = minidom.parse(orngServerFiles.localpath_download("PPI", "allppis.xml"))
        interactions = document.getElementsByTagName("interaction")
        self.interactions = [process(interaction) for interaction in interactions] 
    
    def __iter__(self):
        return iter(self.interactions)
    
    @classmethod
    def download(cls):
        import urllib2, shutil
        src = urllib2.urlopen("http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz")
        dest = orngServerFiles.localpath("PPI", "mppi.gz")
        shutil.copyfileobj(src, open(dest, "wb"))
        
from obiTaxonomy import pickled_cache
@pickled_cache(None, [("PPI", "allppis.xml")], version=1)
def _mips():
    return MIPS()
        

def mips_interactions():
    return list(_mips())