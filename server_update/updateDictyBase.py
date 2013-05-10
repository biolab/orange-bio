##interval:7
from common import *

import sys, os
from gzip import GzipFile
import tempfile
from Orange.bio.obiDicty import DictyBase
import Orange.bio.obiDictyMutants as DictyMutants

tmpdir = tempfile.mkdtemp("dictybase")
base = DictyBase.pickle_data()
filename = os.path.join(tmpdir, "tf")

f = open(filename, 'wb')
f.write(base)
f.close()

dom = DictyBase.domain
fn = DictyBase.filename

try:
    sf_server.create_domain('dictybase')
except:
    pass

print filename

sf_server.upload(dom, fn, filename, title="dictyBase gene aliases",
    tags=DictyBase.tags)
sf_server.unprotect(dom, fn)

shutil.rmtree(tmpdir)


"""
Orange server upload for DictyMutants 
"""

tmpdir_mutants = tempfile.mkdtemp("dictymutants")
base_mutants = DictyMutants.download_mutants()
file_mutants = os.path.join(tmpdir_mutants, "tempMut")

fm = open(file_mutants, "wb")
fm.write(base_mutants)
fm.close()

fm_dom = DictyMutants.domain
fm_name = DictyMutants.pickle_file

print file_mutants

sf_server.upload(fm_dom, fm_name, file_mutants, title="dictyBase mutant collection",
    tags=DictyMutants.tags)
sf_server.unprotect(fm_dom, fm_name)

shutil.rmtree(tmpdir_mutants)

"""
Orange server upload for Dicty mutant gene sets
"""
from Orange.bio.obiGeneSets import dictyMutantSets, update_server_list
import cPickle as pickle

mutant_sets = dictyMutantSets()

tmpdir_sets = tempfile.mkdtemp("dictygenesets")
file_mutants = os.path.join(tmpdir_sets, "tempSets")

with open(file_mutants, "wb") as f:
    pickle.dump(mutant_sets, f, -1)

fm_dom = "gene_sets"
fm_name = "DictyMutant phenotypes"
set_tags = ["Mutant", "phenotype", "dicty"]

print file_mutants

sf_server.upload(fm_dom, fm_name, file_mutants, title="Dictyostelium discoideum mutant phenotypes",
    tags=set_tags)
sf_server.unprotect(fm_dom, fm_name)
update_server_list(sf_server)

shutil.rmtree(tmpdir_sets)

