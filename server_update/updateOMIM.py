##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si


from common import *
from orangecontrib.bio import obiOMIM

import os, sys

path = os.path.join(environ.buffer_dir, "tmp_OMIM")

try:
    os.mkdir(path)
except OSError:
    pass
filename = os.path.join(path, "morbidmap")
obiOMIM.OMIM.download_from_NCBI(filename)

sf_server.upload("OMIM", "morbidmap", filename, title="Online Mendelian Inheritance in Man (OMIM)",
                   tags=["genes", "diseases", "human", "OMIM" "#version:%i" % obiOMIM.OMIM.VERSION])
sf_server.unprotect("OMIM", "morbidmap")


"""
Orange server upload for OMIM morbidmap gene sets
"""
from orangecontrib.bio.geneset import omimGeneSets, register

omim_sets_split = omimGeneSets().split_by_hierarchy()
for omim_sets in omim_sets_split:
    register(omim_sets, sf_server)

