##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from Orange.bio import obiPPI
import urllib2, gzip

from common import *
import re


def get_version():
    from urllib2 import build_opener
    html = build_opener().open('http://www.string-db.org/newstring_cgi/show_download_page.pl').read().decode()
    ver = re.findall("protein\.links\.(v.*?)\.txt\.gz", html, re.DOTALL)[0]
    return ver

version = get_version()
version_id = "#dbversion:%s" % version

force = False # force update

for cl,desc,sfn in [ (obiPPI.STRING, 
                    "STRING Protein interactions (Creative Commons Attribution 3.0 License)", 
                    obiPPI.STRING.FILENAME),
                    (obiPPI.STRINGDetailed, 
                    "STRING Protein interactions (Creative Commons Attribution-Noncommercial-Share Alike 3.0 License)", 
                    obiPPI.STRINGDetailed.FILENAME_DETAILED) ]:

    print cl
    print "current info", sf_server.info("PPI", sfn)

    if force or version_id not in sf_server.info("PPI", sfn)["tags"]:

        filename = sf_local.localpath("PPI",  sfn)

        if os.path.exists(filename):
            os.remove(filename)

        cl.download_data(version)

        gzfile = gzip.GzipFile(filename + ".gz", "wb")
        shutil.copyfileobj(open(filename, "rb"), gzfile)

        sf_server.upload("PPI", sfn, filename + ".gz", 
                       desc,
                       tags=["protein interaction", "STRING", 
                             "#compression:gz", "#version:%s" % cl.VERSION, version_id]
                       )
        sf_server.unprotect("PPI", sfn)
