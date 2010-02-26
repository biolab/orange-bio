from __future__ import division
import urllib
import re
import pylab
import random
import os
import math
import locale
import statc
import numpy as np

import obiTaxonomy
import obiGO as go
import obiProb as op
import orngServerFiles as osf
import obiGene as ge
import obiKEGG as kg


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
[preIDs, premiRNA_lib, preACCtoID, clusters] = __build_lib(premirnafile,0,0,1,1)

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
            objectID = re.sub('miR','mir',objectID)
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
    if accession in preACCtoID:
        return preACCtoID[accession]
    else:
        print "Accession not found."
        return False
    

def __reverseDict(old_dict):
    """
    switch dictionary: keys <--> values; 
    """
    new_dict = {}
    for k,valuesList in old_dict.items():
        for val in valuesList:
            if val != [] and not(val in new_dict):
                new_dict[val]=[]
            if not(k in new_dict[val]):
                new_dict[val].append(k)
    return new_dict


def get_geneMirnaLib(org=None):
    """
    build dictionary gene:[miRNAs]
    """
    mirnaGenes={}
    for m in ids(org):
        mirnaGenes[m] = get_info(m).targets.split(',')
        
    return __reverseDict(mirnaGenes)


def get_GO(mirna_list, annotations, enrichment=False, pval=0.1, goSwitch=True):
    """
    get_GO() takes as input a list of miRNAs of the organism for which the annotations are defined.
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
            mirAnnotations[m] = list(set([an.GO_ID for an in mirna_ann]))
        elif enrichment==True:
            if len(genes):
                resP = annotations.GetEnrichedTerms(genes,aspect='P')
                resC = annotations.GetEnrichedTerms(genes,aspect='C')
                resF = annotations.GetEnrichedTerms(genes,aspect='F')
                res = {}
                res.update(resP)
                res.update(resC)
                res.update(resF)
                tups = [(pVal,go_id) for go_id, (ge,pVal,ref) in res.items()]
                tups.sort()            
                p_correct = op.FDR([p for p,go_id in tups])            
                mirAnnotations[m] = [tups[i][1] for i, p in enumerate(p_correct) if p < pval]
            else:
                mirAnnotations[m]=[]

    if goSwitch:
        return __reverseDict(mirAnnotations)
    else:
        return mirAnnotations
 


def filter_GO(mirna_goid, annotations, treshold=[80,85], reverse=True):    
    """
    removeStopWords() takes as input a dictionary like {mirna:[list of GO_IDs]} and
    remove the most common GO IDs in each list using the TF-IDF criterion.
    Treshold is chosen by making the average of the percentiles introduced.
    """
        
    #mirna_goid = dict(filter(lambda x: x[1], mirna_goid.items()))        
    uniqGO = list(set(reduce(lambda x,y: x+y, mirna_goid.values())))        
    
    goIDF = {}    
    for n,go in enumerate(uniqGO):        
        goIDF[go] = np.log(len(mirna_goid)/len(filter(lambda x: go in x, mirna_goid.values())))        

    
    new_dict={}
    for m,goList in mirna_goid.items():
        TF_IDF ={}
        genes = get_info(m).targets.split(',')
        mirnaAnnotations = [annotations.GetAnnotatedTerms(g).keys() for g in genes]
        for go in goList:
            TF = len(filter(lambda x: go in x, mirnaAnnotations))/len(genes)
            TF_IDF[go] = TF*goIDF[go]
        
        if treshold:
            data = sorted(TF_IDF.values())    
            tresholds = filter(lambda x: x[1]>treshold[0] and x[1]<treshold[1], [(d,statc.percentileofscore(data, d)) for d in list(set(data))])
            t = np.mean([t[0] for t in tresholds])
        else:
            t=0 
           
        new_dict[m] = filter(lambda x: TF_IDF[x] > t, goList) 
    
    if reverse:
        return __reverseDict(new_dict)
    else:
        return new_dict



def get_pathways(mirna_list, organism='hsa', enrichment=False, pVal=0.1, pathSwitch=True):
    """
    get_pathways() takes as input a list of miRNAs and returns a dictionary that has miRNAs as keys
    and pathways IDs as values; if the switch is set on True, 
    it returns a dictionary with pathways IDs as keys and miRNAs as values.
    """
    gmkegg = ge.GMKEGG(organism)
    org = kg.KEGGOrganism(organism)
    gmkegg.set_targets(org.get_genes())     
    
    genes = list(set(reduce(lambda x,y: x+y,[get_info(m).targets.split(',') for m in mirna_list])))
    keggNames = dict([(g,gmkegg.umatch(g)) for g in genes if gmkegg.umatch(g)])
        
    mirnaPathways = {}
    for m in mirna_list:
        kegg_genes = [keggNames[g] for g in get_info(m).targets.split(',') if g in keggNames]
        if enrichment:
            mirnaPathways[m] = [path_id for path_id,(geneList,p,geneNum) in org.get_enriched_pathways_by_genes(kegg_genes).items() if p < pVal]
        else:
            paths = filter(None,[list(org.get_pathways_by_genes([k])) for k in kegg_genes])                   
            if paths:
                mirnaPathways[m] = list(set(reduce(lambda x,y: x+y,paths)))
            else:
                mirnaPathways[m] = []
    
    if pathSwitch:
        return __reverseDict(mirnaPathways)
    else:
        return mirnaPathways
        

    
#######################











      
