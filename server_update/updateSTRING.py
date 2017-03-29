""" update string database """
import re
import gzip

from server_update import *
from server_update.tests.test_STRING import StringTest
from orangecontrib.bio import ppi, taxonomy
from urllib.request import build_opener

DOMAIN = 'PPI'

domain_path = sf_local.localpath(DOMAIN)
downloads = sf_local.localpath('downloaded_files')
temp_path = os.path.join(domain_path, 'temp')

create_folder(domain_path)
create_folder(temp_path)
create_folder(downloads)


def get_version():
    html = build_opener().open('http://string.embl.de/cgi/download.pl').read().decode()
    ver = re.findall("protein\.links\.(v.*?)\.txt\.gz", html, re.DOTALL)[0]
    return ver

version = get_version()
version_id = 'dbversion:{}'.format(version)
force = False  # force update

taxids = ppi.STRING.common_taxids()
desc = "STRING Protein interactions for {name} (Creative Commons Attribution 3.0 License)"


#Sorry, STRING does not know an organism named 'Mycoplasma pneumoniae M129', ' Candida albicans'.
#Note that the database only contains organisms with full sequenced and published genomes.

exclude = ['272634', '5476']
taxids = [idtax for idtax in taxids if idtax not in exclude]

for taxid in taxids:
    dbfilename = ppi.STRING.default_db_filename(taxid)
    basename = os.path.basename(dbfilename)

    TITLE = desc.format(name=taxonomy.name(taxid))
    TAGS = ["protein interaction", "STRING"]
    VERSION = ppi.STRING.VERSION

    if not force and version_id in sf_server.info("PPI", basename)["tags"]:
        continue

    ppi.STRING.init_db(version, taxid, dbfilename=dbfilename, cache_dir=downloads)

    gzfile = gzip.GzipFile(os.path.join(temp_path, basename), "wb")
    shutil.copyfileobj(open(dbfilename, "rb"), gzfile)
    gzfile.close()

    create_info_file(os.path.join(temp_path, basename), title=TITLE, tags=TAGS, version=VERSION,
                     compression='gz', uncompressed=file_size_bytes(dbfilename), dbversion=version_id)


desc_detailed = "STRING Protein interactions for {name} (Creative Commons Attribution-Noncommercial-Share Alike 3.0 License)"


for taxid in taxids:
    print(taxid)
    dbfilename = sf_local.localpath(
        ppi.STRINGDetailed.DOMAIN,
        ppi.STRINGDetailed.FILENAME_DETAILED.format(taxid=taxid)
    )
    basename = os.path.basename(dbfilename)

    TITLE = desc_detailed.format(name=taxonomy.name(taxid))
    TAGS = ["protein interaction", "STRING"]
    VERSION = ppi.STRING.VERSION

    if not force and version_id in sf_server.info("PPI", basename)["tags"]:
        continue

    ppi.STRINGDetailed.init_db(version, taxid, dbfilename=dbfilename, cache_dir=downloads)
    gzfile = gzip.GzipFile(os.path.join(temp_path, basename), "wb")  # gzip the database
    shutil.copyfileobj(open(dbfilename, "rb"), gzfile)
    gzfile.close()

    create_info_file(os.path.join(temp_path, basename), title=TITLE, tags=TAGS, version=VERSION,
                     compression='gz', uncompressed=file_size_bytes(dbfilename), dbversion=version_id)


helper = SyncHelper(DOMAIN, StringTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
