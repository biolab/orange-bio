import gzip
import sqlite3
import time

from server_update import *
from urllib.request import urlopen
from collections import defaultdict
from orangecontrib.bio.gene import homology as obiHomoloGene
from server_update.tests.test_HomoloGene import HomologyTest


DOMAIN = 'HomoloGene'
domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, 'temp')
create_folder(domain_path)
create_folder(temp_path)

file_path = os.path.join(domain_path, 'homologene.data')
compressed_path = os.path.join(temp_path, 'homologene.data')

obiHomoloGene.HomoloGene.download_from_NCBI(file_path)
uncompressed = os.stat(file_path).st_size

with gzip.open(compressed_path, "wb") as f:
    with open(file_path, "rb") as temp_file:
        shutil.copyfileobj(temp_file, f)


TITLE = 'HomoloGene'
TAGS = ["genes", "homologs", "HomoloGene"]
VERSION = obiHomoloGene.HomoloGene.VERSION

print("{} created!".format(compressed_path))
create_info_file(compressed_path, title=TITLE, tags=TAGS, uncompressed=uncompressed,
                 compression='gz', version=VERSION)

# InParanioid Orthologs update

organisms = {"3702": "A.thaliana",
             "9913": "B.taurus",
             "6239": "C.elegans",
             "3055": "C.reinhardtii",
             "7955": "D.rerio",
             "352472": "D.discoideum",
             "7227":  "D.melanogaster",
             "562":  "E.coli",
             # "11103", # Hepatitis C virus
             "9606": "H.sapiens",
             "10090": "M.musculus",
             # "2104",  # Mycoplasma pneumoniae
             "4530": "O.sativa",
             "5833": "P.falciparum",
             # "4754",  # Pneumocystis carinii
             "10116": "R.norvegicus",
             "4932": "S.cerevisiae",
             "4896":  "S.pombe",
             "31033": "T.rubripes"
             # "8355",  # Xenopus laevis
             # "4577",  # Zea mays
             }


def gen(i=0):
    while True:
        yield str(i)
        i += 1


combined_orthologs = []
unique_cluster_id = defaultdict(gen().__next__)
organisms = sorted(organisms.values())

for i, org1 in enumerate(organisms):
    print(i, org1)
    for org2 in organisms[i+1:]:
        filename = "http://inparanoid.sbc.su.se/download/current/Orthologs_OrthoXML/%s/%s-%s.orthoXML" % (org1, org1,
                                                                                                          org2)
        print(filename)
        stream = urlopen(filename)
        orthologs = obiHomoloGene._parseOrthoXML(stream)
        orthologs = [(unique_cluster_id[org1, org2, clid], taxid, gene_symbol)
                     for (clid, taxid , gene_symbol) in orthologs]
        combined_orthologs.extend(orthologs)
        time.sleep(1)

file_path = os.path.join(domain_path, 'InParanoid.sqlite')
compressed_path = os.path.join(temp_path, 'InParanoid.sqlite')

con = sqlite3.connect(file_path)
con.execute("drop table if exists homologs")
con.execute("create table homologs (groupid text, taxid text, geneid text)")
con.execute("create index group_index on homologs(groupid)")
con.execute("create index geneid_index on homologs(geneid)")
con.executemany("insert into homologs values (?, ?, ?)", combined_orthologs)
con.commit()


with open(file_path, 'rb') as file:
    gzfile = gzip.GzipFile(compressed_path, "wb")
    shutil.copyfileobj(file, gzfile)
    gzfile.close()

UNCOMPRESSED = os.stat(file_path).st_size
TITLE = 'InParanoid: Eukaryotic Ortholog Groups'
TAGS = ["genes", "homologs", "orthologs", "InParanoid"]
VERSION = obiHomoloGene.InParanoid.VERSION
create_info_file(compressed_path, title=TITLE, tags=TAGS, version=VERSION,
                 uncompressed=UNCOMPRESSED)

# sync files with remote server
helper = SyncHelper(DOMAIN, HomologyTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
