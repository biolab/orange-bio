""" Taxonomy update """
import tarfile
import bz2


from server_update import *
from orangecontrib.bio import taxonomy
from orangecontrib.bio.ncbi import taxonomy as ncbi_taxonomy
from server_update.tests.test_Taxonomy import TaxonomyTest


DOMAIN = taxonomy.Taxonomy.DOMAIN
FILENAME = taxonomy.Taxonomy.FILENAME
TITLE = "NCBI Taxonomy"
TAGS = ["NCBI", "taxonomy", "organism names", "essential"]

domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, 'temp')
create_folder(domain_path)
create_folder(temp_path)

taxdump_filename = os.path.join(domain_path, "taxdump.tar.gz")
db_filename = os.path.join(domain_path, FILENAME)

ncbi_taxonomy.Taxonomy.download(domain_path)
ncbi_taxonomy.Taxonomy.init_db(db_filename, tarfile.open(taxdump_filename))
create_info_file(db_filename)  # to run tests, we need .info file -> bio.taxonomy.pickled_cache

db_size = os.stat(db_filename).st_size  # store uncompressed database size

with bz2.BZ2File(os.path.join(temp_path, FILENAME), mode="w", compresslevel=9) as f:
    shutil.copyfileobj(open(db_filename, "rb"), f)
create_info_file(os.path.join(temp_path, FILENAME), title=TITLE, tags=TAGS, uncompressed=db_size, compression='bz2')

# sync files with remote server
helper = SyncHelper(DOMAIN, TaxonomyTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
