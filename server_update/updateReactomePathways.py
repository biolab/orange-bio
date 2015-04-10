##interval:7
from common import *

"""
Orange server upload for Cytoband gene sets
"""
from orangecontrib.bio.geneset import reactomePathwaysGeneSets, register

reactome_sets_split = reactomePathwaysGeneSets().split_by_hierarchy()
for pathway_sets in reactome_sets_split:
    register(pathway_sets, sf_server)
