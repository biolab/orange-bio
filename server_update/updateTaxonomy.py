##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from common import *

from Orange.bio import obiTaxonomy

import tarfile

path = os.path.join(environ.buffer_dir, "tmp_Taxonomy")

try:
    os.makedirs(path)
except OSError:
    pass

uncompressedSize = lambda filename: sum(info.size for info in tarfile.open(filename).getmembers())


obiTaxonomy.Taxonomy.ParseTaxdumpFile(outputdir=path)
tFile = tarfile.open(os.path.join(path, "ncbi_taxonomy.tar.gz"), "w:gz")
tFile.add(os.path.join(path, "ncbi_taxonomy.db"), "ncbi_taxonomy.db")
tFile.add(os.path.join(path, "ncbi_taxonomy_inf.db"), "ncbi_taxonomy_inf.db")
tFile.close()


sf_server.upload(
    "Taxonomy", "ncbi_taxonomy.tar.gz",
    os.path.join(path, "ncbi_taxonomy.tar.gz"),
    title="NCBI Taxonomy",
    tags=["NCBI", "taxonomy", "organism names", "essential",
          "#uncompressed:%i" % uncompressedSize(os.path.join(path, "ncbi_taxonomy.tar.gz"))]
)
sf_server.unprotect("Taxonomy", "ncbi_taxonomy.tar.gz")
