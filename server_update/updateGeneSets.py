""" update gene_sets """
import server_update
import orangecontrib.bio.kegg.caching as keggcache

from server_update.tests.test_GeneSets import GeneSetsTest
from orangecontrib.bio.geneset import upload_genesets


keggcache.clear_cache()
upload_genesets(server_update)

helper = server_update.SyncHelper('gene_sets', GeneSetsTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()

