##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from Orange.bio import obiPPI
import urllib2, gzip

from common import *

filename = sf_local.localpath("PPI", obiPPI.STRING.FILENAME)

if False:
    if os.path.exists(filename):
        os.remove(filename)

    obiPPI.STRING.download_data("v9.0")

    gzfile = gzip.GzipFile(filename + ".gz", "wb")
    shutil.copyfileobj(open(filename, "rb"), gzfile)

    sf_server.upload("PPI", obiPPI.STRING.FILENAME, filename + ".gz", 
                       "STRING Protein interactions (Creative Commons Attribution 3.0 License)",
                       tags=["protein interaction", "STRING", 
                             "#compression:gz", "#version:%s" % obiPPI.STRING.VERSION]
                       )
    sf_server.unprotect("PPI", obiPPI.STRING.FILENAME)

# The second part
filename = sf_local.localpath("PPI", obiPPI.STRINGDetailed.FILENAME_DETAILED)

if os.path.exists(filename):
    os.remove(filename)

obiPPI.STRINGDetailed.download_data("v9.0")

gzfile = gzip.GzipFile(filename + ".gz", "wb")
shutil.copyfileobj(open(filename, "rb"), gzfile)

sf_server.upload("PPI", obiPPI.STRINGDetailed.FILENAME_DETAILED, filename + ".gz", 
                   "STRING Protein interactions (Creative Commons Attribution-Noncommercial-Share Alike 3.0 License)" ,
                   tags=["protein interaction", "STRING",
                         "#compression:gz", "#version:%s" % obiPPI.STRINGDetailed.VERSION]
                   )
sf_server.unprotect("PPI", obiPPI.STRINGDetailed.FILENAME_DETAILED)
    
