##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

from common import *

from Orange.bio import obiTaxonomy
import Orange.utils.serverfiles as orngServerFiles
import Orange.utils.environ

import tarfile
import socket

path = os.path.join(environ.buffer_dir, "tmp_Taxonomy")
u = obiTaxonomy.Update(local_database_path=path)

uncompressedSize = lambda filename: sum(info.size for info in tarfile.open(filename).getmembers())

if u.IsUpdatable(obiTaxonomy.Update.UpdateTaxonomy, ()):
    for i in range(3):
        try:
            u.UpdateTaxonomy()
            break
        except socket.timeout, ex:
            print ex
            pass
    sf_server.upload("Taxonomy", "ncbi_taxonomy.tar.gz", os.path.join(path, "ncbi_taxonomy.tar.gz"), title ="NCBI Taxonomy",
                       tags=["NCBI", "taxonomy", "organism names", "essential", "#uncompressed:%i" % uncompressedSize(os.path.join(path, "ncbi_taxonomy.tar.gz"))])
    sf_server.unprotect("Taxonomy", "ncbi_taxonomy.tar.gz")
