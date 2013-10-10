#!/usr/bin/env python

from collections import defaultdict
from Orange.bio.obiGeneAtlas import run_simple_query
from Orange.bio.obiGeneSets import GeneSets, GeneSet
import time

def display_string(string):
    return string.capitalize().replace("_", " ")

regulation = "updown"
organism   = "Homo sapiens"
condition  = "organism_part"
org_code   = "9606"
max_pvalue = 1e-5

query = run_simple_query(regulation=regulation, organism=organism, condition=condition, start=0, rows=5) # Do a short pre-query to see how many genes there are
total_genes = query["totalResults"]

no_of_genes = 50 # 50 genes per query used, since the HTTP requests don't seem stable with more than 50 genes at a time.
no_of_pages = total_genes / no_of_genes  # Gene Expression Atlas API only permits access to [50,200] genes at a time.

print total_genes
#exit()

i=0
sets = defaultdict(list)

for start in range(no_of_pages+1):
    start = time.time()
    query = run_simple_query(regulation=regulation, organism=organism, condition=condition, start=start*no_of_genes, rows=no_of_genes)
    for result in query["results"]:
        i += 1
        print result["gene"]["name"] + "\t" + str(i) # For printing out gene names and the current number of genes during debugging
        for exp in result["expressions"]: 
            diff_exp = [e for e in exp["experiments"] if e["pvalue"] <= max_pvalue] # Use only genes that are significantly diff. expressed
            if diff_exp:
                try:
                    sets[exp["ef"], exp["efv"]].append(result["gene"]["name"]) # The Gene Expression Atlas entries are not consistent. 
                except:
                    sets[exp["efoTerm"], exp["efoId"]].append(result["gene"]["name"])
    print time.time() - start
        
gene_sets = []
for (ef, efv), genes in sets.items():
    ef_display = display_string(ef)
    gs = GeneSet(genes, "Diff. expressed in %s=%s." % (ef_display, efv), id=ef + ":" + efv,
                 description="Diff. expressed in %s=%s." % (ef_display, efv),
                 link="http://www.ebi.ac.uk/gxa/qrs?specie_0={organism}&gprop_0=&gnot_0=&gval_0=&fact_1=&fexp_1=UPDOWN&fmex_1=&fval_1=%22{efv}%22+&view=hm".format( \
                         organism = "+".join(organism.lower().split()), efv = "+".join(efv.lower().split())), organism=org_code, hierarchy=("Gene Expression Atlas", display_string(condition)))
    gene_sets.append(gs)

final_set = GeneSets(gene_sets)

print final_set
exit()

from Orange.bio.obiGeneSets import register
import Orange.utils.serverfiles as serverfiles
import sys

try:
    sf_server = serverfiles.ServerFiles(sys.argv[1], sys.argv[2])
except:
    print "argv[1] = username, argv[2] = password"
    exit()

set_split = final_set.split_by_hierarchy()

for s in set_split:
    register(s, sf_server)

print "Gene sets successfully registered..."
