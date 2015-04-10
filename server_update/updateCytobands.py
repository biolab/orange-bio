##interval:7
from common import *

"""
Orange server upload for Cytoband gene sets
"""
from orangecontrib.bio.geneset import cytobandGeneSets, register

cytoband_sets_split = cytobandGeneSets().split_by_hierarchy()
for band_sets in cytoband_sets_split:
    register(band_sets, sf_server)
