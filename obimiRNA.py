

import urllib
import re
import pylab
import random
import os
import math
import locale

import obiTaxonomy
import orngServerFiles as osf


mirnafile = osf.localpath_download('miRNA','miRNA.txt')
premirnafile = osf.localpath_download('miRNA','premiRNA.txt')

###################################################################

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
[IDs, LABELS, miRNA_lib, mat_toPre] = __build_lib(mirnafile, 1,1,0,0)
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


#######################

                 