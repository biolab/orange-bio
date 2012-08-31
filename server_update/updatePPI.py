##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from Orange.bio import obiPPI
import Orange.utils.serverfiles as orngServerFiles
import os, sys, shutil, urllib2, tarfile
from getopt import getopt

opt = dict(getopt(sys.argv[1:], "u:p:", ["user=", "password="])[0])

username = opt.get("-u", opt.get("--user", "username"))
password = opt.get("-p", opt.get("--password", "password"))

serverFiles = orngServerFiles.ServerFiles(username, password)

try:
    os.mkdir(orngServerFiles.localpath("PPI"))
except OSError:
    pass

try:
    serverFiles.create_domain("PPI")
except Exception, ex:
    print ex

if True:
    obiPPI.MIPS.download()

    filename = orngServerFiles.localpath("PPI", "mppi.gz")
    serverFiles.upload("PPI", "allppis.xml", filename, "MIPS Protein interactions",
                       tags=["protein interaction", "MIPS", "#compression:gz", "#version:%i" % obiPPI.MIPS.VERSION]
                       )
    serverFiles.unprotect("PPI", "allppis.xml") 

if True:
    obiPPI.BioGRID.download_data("http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.1.91/BIOGRID-ALL-3.1.91.tab2.zip") #replace with the newest version

    sfn = obiPPI.BioGRID.SERVER_FILE

    filename = orngServerFiles.localpath("PPI", sfn)

    import gzip
    gz = gzip.GzipFile(filename + ".gz", "wb")
    gz.write(open(filename, "rb").read())
    gz.close()

    serverFiles.upload("PPI", sfn, filename + ".gz", 
        title="BioGRID Protein interactions", 
        tags=["protein interaction", "BioGrid", "#compression:gz", "#version:%s" % obiPPI.BioGRID.VERSION]
        )
    serverFiles.unprotect("PPI", sfn)

