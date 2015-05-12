##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from common import *

import tarfile
import bz2

from orangecontrib.bio import taxonomy
from orangecontrib.bio.ncbi import taxonomy as ncbi_taxonomy


DOMAIN, FILENAME = taxonomy.Taxonomy.DOMAIN, taxonomy.Taxonomy.FILENAME

path = os.path.join(environ.buffer_dir, "tmp_Taxonomy")

try:
    os.makedirs(path)
except OSError:
    pass

taxdump_filename = os.path.join(path, "taxdump.tar.gz")
db_filename = os.path.join(path, FILENAME)
bz2_filename = os.path.join(path, FILENAME + ".bz2")

ncbi_taxonomy.Taxonomy.download(path)
ncbi_taxonomy.Taxonomy.init_db(db_filename, tarfile.open(taxdump_filename))

db_size = os.stat(db_filename).st_size

with bz2.BZ2File(bz2_filename, mode="w", compresslevel=9) as f:
    shutil.copyfileobj(open(db_filename, "rb"), f)

sf_server.upload(
    DOMAIN, FILENAME,
    bz2_filename,
    title="NCBI Taxonomy",
    tags=["NCBI", "taxonomy", "organism names", "essential",
          "#uncompressed:%i" % db_size, "#compression:bz2"]
)
sf_server.unprotect(DOMAIN, FILENAME)
