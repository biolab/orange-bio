#!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from Orange.bio import ppi, taxonomy
import urllib2, gzip

from common import *
import re

tmp_path = os.path.join(environ.buffer_dir, "tmp_pack_STRING")
try:
    os.makedirs(tmp_path)
except Exception:
    pass

def get_version():
    from urllib2 import build_opener
    html = build_opener().open('http://string-db.org/newstring_cgi/show_download_page.pl').read().decode()
    ver = re.findall("protein\.links\.(v.*?)\.txt\.gz", html, re.DOTALL)[0]
    return ver

version = get_version()
version_id = "#dbversion:%s" % version

force = False # force update

taxids = ppi.STRING.common_taxids()
desc = "STRING Protein interactions for {name} (Creative Commons Attribution 3.0 License)"

for taxid in taxids:
    dbfilename = ppi.STRING.default_db_filename(taxid)
    basename = os.path.basename(dbfilename)

    if not force and version_id in sf_server.info("PPI", basename)["tags"]:
        continue

    ppi.STRING.init_db(version, taxid, dbfilename=dbfilename, cache_dir=tmp_path)
    gzfile = gzip.GzipFile(dbfilename + ".gz", "wb")  # gzip the database
    shutil.copyfileobj(open(dbfilename, "rb"), gzfile)
    gzfile.close()
    sf_server.remove("PPI", basename)
    sf_server.upload(
        "PPI", basename, dbfilename + ".gz",
        desc.format(name=taxonomy.name(taxid)),
        tags=["protein interaction", "STRING",
              "#compression:gz",
              "#version:%s" % ppi.STRING.VERSION,
              version_id]
    )
    sf_server.unprotect("PPI", basename)

desc_detailed = "STRING Protein interactions for {name} (Creative Commons Attribution-Noncommercial-Share Alike 3.0 License)"

force = False  # force updatea

for taxid in taxids:
    dbfilename = sf_local.localpath(
        ppi.STRINGDetailed.DOMAIN,
        ppi.STRINGDetailed.FILENAME_DETAILED.format(taxid=taxid)
    )
    basename = os.path.basename(dbfilename)

    if not force and version_id in sf_server.info("PPI", basename)["tags"]:
        continue

    ppi.STRINGDetailed.init_db(version, taxid, dbfilename=dbfilename, cache_dir=tmp_path)
    gzfile = gzip.GzipFile(dbfilename + ".gz", "wb")  # gzip the database
    shutil.copyfileobj(open(dbfilename, "rb"), gzfile)
    gzfile.close()
    sf_server.remove("PPI", basename)
    sf_server.upload(
        "PPI", basename, dbfilename + ".gz",
        desc_detailed.format(name=taxonomy.name(taxid)),
        tags=["protein interaction", "STRING",
              "#compression:gz",
              "#version:%s" % ppi.STRING.VERSION,
              version_id]
    )
    sf_server.unprotect("PPI", basename)
