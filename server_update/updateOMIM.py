##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from common import *
from Orange.bio import obiOMIM

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
