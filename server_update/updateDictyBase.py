""" Dictybase update """
from server_update import *
from server_update.tests.test_Dicty import DictyBaseTest, DictyMutantsTest
from orangecontrib.bio.dicty import phenotypes as DictyMutants
from orangecontrib.bio.dicty import DictyBase

"""
Update dictyBase
"""
DOMAIN = DictyBase.domain
FILENAME = DictyBase.filename
TITLE = 'dictyBase gene aliases'
TAGS = DictyBase.tags
base = DictyBase.pickle_data()
create_folder(sf_local.localpath(DOMAIN))
localfile = sf_local.localpath(DOMAIN, FILENAME)

print("updating dictBase ")
with open(localfile, 'wb') as f:
    f.write(base)
    f.close()

create_info_file(localfile, title=TITLE, tags=TAGS)

helper = SyncHelper(DOMAIN, DictyBaseTest)
helper.run_tests()
helper.sync_files()

print("DictyBase updated\n")


"""
Update DictyMutants
"""
base_mutants = DictyMutants.download_mutants()
DOMAIN = DictyMutants.domain
FILENAME = DictyMutants.pickle_file
TITLE = 'dictyBase mutant phenotypes'
TAGS = DictyMutants.tags
create_folder(sf_local.localpath(DOMAIN))
localfile = sf_local.localpath(DOMAIN, FILENAME)


print("updating dictBase mutants")
with open(localfile, 'wb') as f:
    f.write(base_mutants)
    f.close()

create_info_file(localfile, title=TITLE, tags=TAGS)

helper = SyncHelper(DOMAIN, DictyMutantsTest)
helper.run_tests()
helper.sync_files()
print("DictyBase mutant updated")

helper.remove_update_folder()
