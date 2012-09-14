##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from Orange.bio import obiPPI
import urllib2, tarfile

from common import *

try:
    os.mkdir(sf_local.localpath("PPI"))
except OSError:
    pass

try:
    sf_server.create_domain("PPI")
except Exception, ex:
    print ex

if True:
    obiPPI.MIPS.download()

    filename = sf_local.localpath("PPI", "mppi.gz")
    sf_server.upload("PPI", "allppis.xml", filename, "MIPS Protein interactions",
                       tags=["protein interaction", "MIPS", "#compression:gz", "#version:%i" % obiPPI.MIPS.VERSION]
                       )
    sf_server.unprotect("PPI", "allppis.xml") 

if True:
    obiPPI.BioGRID.download_data("http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.1.91/BIOGRID-ALL-3.1.91.tab2.zip") #replace with the newest version

    sfn = obiPPI.BioGRID.SERVER_FILE

    filename = sf_local.localpath("PPI", sfn)

    import gzip
    gz = gzip.GzipFile(filename + ".gz", "wb")
    gz.write(open(filename, "rb").read())
    gz.close()

    sf_server.upload("PPI", sfn, filename + ".gz", 
        title="BioGRID Protein interactions", 
        tags=["protein interaction", "BioGrid", "#compression:gz", "#version:%s" % obiPPI.BioGRID.VERSION]
        )
    sf_server.unprotect("PPI", sfn)

