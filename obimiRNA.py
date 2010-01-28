
import urllib
import re
import pylab
import random
import os
import math
import locale

import obiTaxonomy
import obiGO as go
import obiProb as op
import orngServerFiles as osf


mirnafile = osf.localpath_download('miRNA','miRNA.txt')
premirnafile = osf.localpath_download('miRNA','premiRNA.txt')

################################################################################################################
################################################################################################################

def __build_lib(filename, labels=True, MATtoPRE=True, ACCtoID=True, clust=False):
    """
    build_lib() function takes as input a filename
    and gives as output some variables there will be used in
    the module.
    """
    content = [l.rstrip() for l in open(filename).readlines()][1:]
    to_return = []
    
    ids = [l.split('\t')[0] for l in content]
    to_return.append(ids)
    
    if labels: 
        to_return.append(list(set(i.split('-')[0] for i in ids)))
    
    elements = [line.rstrip().split('\t') for line in content] 
    to_return.append(dict((elem[0],elem[1:]) for elem in elements))
    
    if MATtoPRE:
        to_return.append(dict([(e[0],e[3]) for e in elements]))
    
    if ACCtoID:
        to_return.append(dict([(e[1],e[0]) for e in elements]))
        
    if clust:
        to_return.append(dict([(e[0],e[5]) for e in elements]))
        
    
    return to_return

### library:    
[IDs, LABELS, miRNA_lib, mat_toPre, ACCtoID] = __build_lib(mirnafile, 1,1,1,0)
[preIDs, premiRNA_lib,clusters] = __build_lib(premirnafile,0,0,0,1)

fromTaxo = {3702:'ath', 9913:'bta', 6239:'cel', 3055:'cre', 7955:'dre',\
             352472:'ddi', 7227:'dme', 9606:'hsa', 10090:'mmu', 4530:'osa',\
              10116:'rno', 8355:'xla', 4577:'zma'}

toTaxo = dict(zip(fromTaxo.values(),fromTaxo.keys()))

num_toClusters = {}
clusters_toNum = {}
n=0
for k,v in clusters.items():
    if v !='None':
        g = v.split(',')
        g.append(k)        
        group = sorted(g)
        if not(group in num_toClusters.values()):
             num_toClusters[n] = group
             for e in group:
                 clusters_toNum[e]=n
             n += 1



class miRNAException(Exception):
    pass

def ids(taxid=None):
    """
    ids() functions takes an organism identifier (for human can be 9606 or hsa)
    and returns all the miRNAs related to that organism. If there 
    is no argument, it returns all the miRNAs in the library.
    """
    
    if not(taxid):
        return IDs
    else:
        taxid = fromTaxo.get(taxid, taxid)
        
        if (taxid in toTaxo):
            return [e for e in IDs if e.split('-')[0]==taxid]
                
        else:
            raise miRNAException("ids() Error: Check the input value.")

            
        
class mat_miRNA:
    pass

class pre_miRNA:
    pass
 
      
def get_info(objectID,type='mat'):
        """
        get_info() function takes a miRNA identifier as input
        and returns a miRNA object.
        """
        if type == 'mat':
            objectID = re.sub('mir','miR',objectID)
            if objectID in IDs:
                attr = [line.rstrip() for line in open(mirnafile).readlines()][0].split('\t')
            
                to_return = mat_miRNA()
                setattr(to_return, attr[0], objectID)
            
                for n,a in enumerate(attr[1:]):
                    setattr(to_return, a, miRNA_lib[objectID][n])
            
                return to_return
            else:
                raise miRNAException("get_info() Error: Check the input value.")
            
        elif type == 'pre':
            if objectID in preIDs:
                attr = [line.rstrip() for line in open(premirnafile).readlines()][0].split('\t')
            
                to_return = pre_miRNA()
                setattr(to_return, attr[0], objectID)
            
                for n,a in enumerate(attr[1:]):
                    setattr(to_return, a, premiRNA_lib[objectID][n])
                            
                return to_return
            else:
                raise miRNAException("get_info() Error: Check the input value.")
        else:
           raise miRNAException("get_info() Error: Check type value.")

                    
def cluster(clusterID, type='name'):
    """
    cluster() function take a cluster identifier or a premiRNA
    and return the list of premiRNAs clustered together."
    """
    if type=='name':
        if clusterID in clusters:
            return clusters[clusterID]
        else:
            raise miRNAException("cluster() Error: ClusterID not found in premiRNA names.")
    
    elif type=='num':
        if clusterID in num_toClusters:
            return num_toClusters[clusterID]
        else:
            raise miRNAException("cluster() Error: ClusterID not found in clusters' list.")
    else:
        raise miRNAException("cluster() Error: Check the input value.")


def fromACC_toID(accession):
    """
    fromACC_toID() takes a miRNA accession number
    and returns a miRNA id.
    """
    if accession in ACCtoID:
        return ACCtoID[accession]
    else:
        print "Accession not found."
        return False
    
    
def __getGO_dict(mirna_dict):
    """
    switch miRNA --> GO_IDs dictionary to GO_ID --> miRNAs. 
    """
    
    GO_mirna_lib = {}
    for m,GO_IDs in mirna_dict.items():
        for go_id in GO_IDs:
            if go_id != [] and not(go_id in GO_mirna_lib):
                GO_mirna_lib[go_id]=[]
            GO_mirna_lib[go_id].append(m)
    
    return GO_mirna_lib


def get_GO(mirna_list, annotations, enrichment=False, pval=0.1, goSwitch=True):
    """
    get_GO() takes as input a list of miRNAs only for the organism for which the annotations are defined.
    If goSwitch is False, get_GO() returns a dictionary that has miRNAs as keys and GO IDs as values;
    in the other case it returns a dictionary with GO IDs as keys and miRNAs as values.
    """
    
    import obiGene
    genematcher = obiGene.matcher([obiGene.GMGO(annotations.taxid)] + \
        ([obiGene.GMDicty()] if annotations.taxid == "352472"  else []))
    genematcher.set_targets(annotations.geneNames)
    
    mirna_list = list(set(mirna_list))
    mirAnnotations = {}
    
    for m in mirna_list:
        genes = get_info(m).targets.split(',')
        genes = filter(None, map(genematcher.umatch, genes))
        
        if enrichment==False:
            mirna_ann=[]
            for gene in genes:
                gene_annotations = annotations.geneAnnotations[gene]
                for ga in gene_annotations:
                    mirna_ann.append(ga)
            mirAnnotations[m] = [an.GO_ID for an in list(set(mirna_ann))]
        elif enrichment==True:
            res = annotations.GetEnrichedTerms(genes)
            tups = [(pVal,go_id) for go_id, (ge,pVal,ref) in res.items()]
            tups.sort()            
            p_correct = op.FDR([p for p, go_id in tups])            
            mirAnnotations[m] = [tups[i][1] for i, p in enumerate(p_correct) if p < pval]

    if goSwitch:
        return __getGO_dict(mirAnnotations)
    else:
        return mirAnnotations
 
    
#######################











      
