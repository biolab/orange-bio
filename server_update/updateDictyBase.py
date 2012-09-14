##interval:7
from common import *

import sys, os
from gzip import GzipFile
import tempfile
from Orange.bio.obiDicty import DictyBase

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
