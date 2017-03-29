import os

from server_update import create_info_file, sf_local, create_folder, SyncHelper
from server_update.tests.test_OMIM import OMIMTest
from orangecontrib.bio import omim


DOMAIN = 'OMIM'
domain_path = sf_local.localpath(DOMAIN)
create_folder(domain_path)


FILENAME = 'morbidmap'
TITLE = 'Online Mendelian Inheritance in Man (OMIM)'
TAGS = ['genes', 'diseases', 'human', 'OMIM']
VERSION = omim.OMIM.VERSION


file_path = os.path.join(domain_path, "morbidmap")
print('downloading file ...')
omim.OMIM.download_from_NCBI(file_path)
print('file created ...')
create_info_file(file_path, title=TITLE, tags=TAGS, version=VERSION)

helper = SyncHelper(DOMAIN, OMIMTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
