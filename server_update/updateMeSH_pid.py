from common import *

from collections import defaultdict

import os
import os.path

import orangecontrib.bio.geneset
orangecontrib.bio.geneset.update_server_list(sf_server)

chemical_term = defaultdict(set)

path = os.path.join(environ.buffer_dir, "tmp_mesh")
try:
    os.makedirs(path)
except Exception:
    pass

def download_file(url, local_filename):
    from urllib import urlopen
    urlopen(url)
    with open(local_filename, 'wb') as f:
        shutil.copyfileobj(urlopen(url), f)
    return local_filename

year = 2015

print path
download_file("ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/d%d.bin" % year, os.path.join(path, "d%d.bin" % year))
download_file("ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/c%d.bin" % year, os.path.join(path, "c%d.bin" % year))
download_file("ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-MeSH", os.path.join(path, "CID-MeSH.txt")) 

f = open(os.path.join(path, "c2015.bin"), "rt")
cur = None
for l in f:
    l = l.strip()
    if l == "*NEWRECORD":
        cur = {}
    if l.startswith("NM = "):
        cur["NM"] = l[5:]
    elif l.startswith("HM = ") or l.startswith("PA = "):
        chemical_term[cur["NM"]].add(l[5:])

term_tree = defaultdict(set)
term_ms = defaultdict(str)
tree_term = {}
term_pa = defaultdict(set)

f = open(os.path.join(path, "d2015.bin"), "rt")
cur = None
for l in f:
    l = l.strip()
    if l == "*NEWRECORD":
        cur = {}
    if l.startswith("MH = "):
        cur["MH"] = l[5:]
    elif l.startswith("PA = "):
        term_pa[cur["MH"]].add(l[5:])
    elif l.startswith("MN = "):
        term_tree[cur["MH"]].add(l[5:])
        if l[5:] in tree_term:
            errorr
        else:
            tree_term[l[5:]] = cur["MH"]
    elif l.startswith("MS = "):
        if cur["MH"] in term_ms:
            errorr
        else:
            term_ms[cur["MH"]] = l[5:]

print term_tree.items()[:5]

term_pids = defaultdict(set)

def addterm(term, pid):
    endpoints = term_tree[term]
    for endp in endpoints:
        if endp.startswith("D"): #chemicals
            parts = endp.split(".")
            for ai in range(1, len(parts) + 1):
                treeid = ".".join(parts[:ai])
                term_pids[tree_term[treeid]].add(pid)
                pas = term_pa[tree_term[treeid]]
                for pa in pas:
                    addterm(pa, pid)
                    


f = open(os.path.join(path, "CID-MeSH.txt"), "rt")
for l in f:
    l = l.strip()
    t = l.split('\t')
    pid, chem = t[0], t[1]
    terms = chemical_term.get(chem, set([chem]))
    for term in terms:
        addterm(term, pid)

for a,b in sorted(term_pids.items(), key=lambda x: len(x[1])):
    print term_tree[a], a, len(b)

from orangecontrib.bio.geneset import GeneSet, GeneSets

gss = set()
for a,b in term_pids.items():
    id_ = sorted(term_tree[a])[0]
    g = GeneSet(genes=b, name=a, id=id_, \
        description=term_ms[a], link="http://www.ncbi.nlm.nih.gov/mesh/" + id_, organism=None, \
        hierarchy=("MeSH", "Chemicals"))
    if not a.startswith("Supplement"):
        gss.add(g)

gss = GeneSets(gss)

orangecontrib.bio.geneset.register(gss, sf_server)

#gssgmt = "chem_mesh2.gmt"
#with open(gssgmt, "wt") as f:
#    for g in gss:
#        f.write('%s\t%s\t%s\n' % (g.name, g.description + 
#            " [http://www.ncbi.nlm.nih.gov/mesh/" + g.id  + "]" , '\t'.join(g.genes)))
                    
    
