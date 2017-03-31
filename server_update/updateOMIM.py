import os

from server_update import *
from server_update.tests.test_OMIM import OMIMTest
from orangecontrib.bio import omim


DOMAIN = 'OMIM'
FILENAME = 'morbidmap'
TITLE = 'Online Mendelian Inheritance in Man (OMIM)'
TAGS = ['genes', 'diseases', 'human', 'OMIM']
VERSION = omim.OMIM.VERSION
helper = SyncHelper(DOMAIN, OMIMTest)

if not http_last_modified('http://ftp.ncbi.nih.gov/repository/OMIM/ARCHIVE/morbidmap') > sf_last_modified(DOMAIN, FILENAME):
    helper.remove_update_folder()
    sys.exit(up_to_date)


domain_path = sf_local.localpath(DOMAIN)
create_folder(domain_path)

file_path = os.path.join(domain_path, "morbidmap")
print('downloading file ...')
omim.OMIM.download_from_NCBI(file_path)
print('file created ...')
create_info_file(file_path, title=TITLE, tags=TAGS, version=VERSION)

helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
