from __future__ import absolute_import, division

from collections import defaultdict
import math, os, random, re, urllib

from Orange.orng import orngServerFiles as osf
import statc

from . import gene as ge, go, kegg as kg, utils, taxonomy as obiTaxonomy

op = utils.stats

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

def open_microcosm(org="mus_musculus", version="v5"):
    """ Open the miRna targets from the EBI microcosm site. 
    """
    import urllib2, zipfile
    import StringIO
    stream = urllib2.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/microcosm/{version}/arch.{version}.txt.{org}.zip".format(org=org, version=version))
    contents = stream.read()
    stream = StringIO.StringIO(contents)
    return stream
    
def parse_targets_microcosm_v5(file, max_pvalue=None, min_score=None):
    """ Parse miRna targets (version 'v5.txt' format file from microcosm) 
    """
    
    ids = set() #mirna ids
    labels = set() # org codes (i.e hsa, mmu ..)
    #["mirna accession", "seq", "mirna id", "targets"]
    mirna_lib = defaultdict(lambda: ["", "", "", set()])
    mat_to_pre = {}
    if isinstance(file, basestring):
        file = open(file, "rb")
    import csv
    reader = csv.reader(file, delimiter="\t")
    for line in reader:
        if not line or line[0].startswith("##"):
            continue
        
        seq = line[1]
        org = seq.split("-", 1)[0] #SEQ
        ids.add(seq)
        labels.add(org)
        
        mirna_lib[seq] #In case some mirna dont have any targets
        
        if max_pvalue is not None:
            if float(line[10]) > max_pvalue:
                continue
        if min_score is not None:
            if float(line[9]) < min_score:
                continue
            
        mirna_lib[seq][-1].add(line[11]) #TRANSCRIPT_ID
        
    for key, value in mirna_lib.iteritems():
        value[-1] = ",".join(sorted(value[-1]) or ["None"])
        
    ids = sorted(ids)
    labels = sorted(labels)
    mirna_lib = dict(mirna_lib)
    return ids, labels, mirna_lib, {}, {}


def load_miRNA_microCosm(org="mus_musculus", max_pvalue=None, min_score=None):
    """ Load miRNA's from microcosm into the global scope (currently
    only Mus musculus is supported)
    
    """
    global IDs, LABELS, miRNA_lib, mat_toPre, ACCtoID
    global preIDs, premiRNA_lib, preACCtoID, clusters
    global num_toClusters,  clusters_toNum
    
    file = osf.localpath_download("miRNA", "v5.txt.{org}".format(org=org))
    [IDs, LABELS, miRNA_lib, mat_toPre, ACCtoID] = parse_targets_microcosm_v5(file,
                                max_pvalue=max_pvalue, min_score=min_score)
    [preIDs, premiRNA_lib, preACCtoID, clusters] = [], {}, {}, {}
    num_toClusters, clusters_toNum = {}, {}

    
load_miRNA = load_miRNA_microCosm


def load_miRNA_TargetScan():
    """ This loads miRNAs from miRBase and targets from TargetScan.
    Will also load pre-miRNAs
    """
    global IDs, LABELS, miRNA_lib, mat_toPre, ACCtoID
    global preIDs, premiRNA_lib, preACCtoID, clusters
    global num_toClusters,  clusters_toNum
    
    [IDs, LABELS, miRNA_lib, mat_toPre, ACCtoID] = __build_lib(mirnafile, 1,1,1,0)
    [preIDs, premiRNA_lib, preACCtoID, clusters] = __build_lib(premirnafile,0,0,1,1)
    
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
     
### library:
load_miRNA()

#[IDs, LABELS, miRNA_lib, mat_toPre, ACCtoID] = __build_lib(mirnafile, 1,1,1,0)
#[preIDs, premiRNA_lib, preACCtoID, clusters] = __build_lib(premirnafile,0,0,1,1)
#
#num_toClusters = {}
#clusters_toNum = {}
#n=0
#for k,v in clusters.items():
#    if v !='None':
#        g = v.split(',')
#        g.append(k)        
#        group = sorted(g)
#        if not(group in num_toClusters.values()):
#             num_toClusters[n] = group
#             for e in group:
#                 clusters_toNum[e]=n
#             n += 1


fromTaxo = {3702:'ath', 9913:'bta', 6239:'cel', 3055:'cre', 7955:'dre',\
             352472:'ddi', 7227:'dme', 9606:'hsa', 10090:'mmu', 4530:'osa',\
              10116:'rno', 8355:'xla', 4577:'zma'}

toTaxo = dict(zip(fromTaxo.values(), fromTaxo.keys()))


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
            from . import obiKEGG
            raise obiTaxonomy.UnknownSpeciesIdentifier(taxid)
        
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
    
    from . import gene as obiGene
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
 


def filter_GO(mirna_goid, annotations, treshold=0.04, reverse=True):    
    """
    filter_GO() takes as input a dictionary like {mirna:[list of GO_IDs]} and
    remove the most common GO IDs in each list using the TF-IDF criterion.
    """       
    uniqGO = list(set(reduce(lambda x,y: x+y, mirna_goid.values())))        
    
    goIDF = {}    
    for n,go in enumerate(uniqGO):        
        goIDF[go] = np.log(len(mirna_goid)/len(filter(lambda x: go in x, mirna_goid.values())))        

    data = []
    new_dict={}
    for m,goList in mirna_goid.items():
        TF_IDF ={}
        genes = get_info(m).targets.split(',')
        mirnaAnnotations = [annotations.GetAnnotatedTerms(g).keys() for g in genes]
        for go in goList:
            TF = len(filter(lambda x: go in x, mirnaAnnotations))/len(genes)
            TF_IDF[go] = TF*goIDF[go]
        if treshold:
            t = treshold
        else:
            t=0
        new_dict[m] = filter(lambda x: TF_IDF[x] > t, goList)   
        #data.append(TF_IDF.values())
       
#    if treshold:
#        data = sorted(reduce(lambda x,y: x+y, data))    
#        tresholds = filter(lambda x: x[1]>treshold[0] and x[1]<treshold[1], [(d,statc.percentileofscore(data, d)) for d in list(set(data))])
#        t = np.mean([t[0] for t in tresholds])
#        print t
#    else:
#        t=0
#    
#    new_dict={}     
#    for m,goList in mirna_goid.items():       
#        new_dict[m] = filter(lambda x: TF_IDF[x] > t, goList) 
    
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
            mirnaPathways[m] = [path_id for path_id,(geneList,p,geneNum) in org.get_enriched_pathways(kegg_genes).items() if p < pVal]
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
        

def removeOldMirnas(mirna_list, getOnlyMature=False):
    """
    removeOldMirnas() takes a list of miRNAs as input and
    divides them in two lists, accordin if they're still present
    on miRBase or not.
    """
    old = []
    for m in mirna_list:
        try:
            mat = get_info(m)
        except Exception:
            if getOnlyMature:
                old.append(m)
            else:                
                try:
                    pre = get_info(m,type='pre')
                except Exception:
                    old.append(m)
    return [old, filter(lambda x: not(x in old), mirna_list)]
    
#######################











      
