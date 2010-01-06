
import urllib, re, pylab, random, os, math, locale, obiTaxonomy
import orngServerFiles as osf


mirnafile = osf.localpath_download('miRNA','miRNA.txt')
premirnafile = osf.localpath_download('miRNA','premiRNA.txt')


IDs = [line.rstrip().split('\t')[0] for line in open(mirnafile).readlines()][1:]

LABELS = list(set(i.split('-')[0] for i in IDs))

miRNA_lib = dict((elements[0],elements[1:]) for elements in [line.rstrip().split('\t') for line in open(mirnafile).readlines()])

preIDs = [line.rstrip().split('\t')[0] for line in open(premirnafile).readlines()][1:]

premiRNA_lib = dict((elements[0],elements[1:]) for elements in [line.rstrip().split('\t') for line in open(premirnafile).readlines()])

fromTaxo = {3702:'ath', 9913:'bta', 6239:'cel', 3055:'cre', 7955:'dre', 352472:'ddi', 7227:'dme', 9606:'hsa', 10090:'mmu', 4530:'osa', 10116:'rno', 8355:'xla', 4577:'zma'}

mat_toPre = dict([(row.rstrip().split('\t')[0],row.rstrip().split('\t')[-1].split(',')) for row in open(mirnafile).readlines()[1:]])

clusters = dict([(c_line.split('\t')[0], c_line.split('\t')[5]) for c_line in open(premirnafile).readlines()[1:]])

num_clusters = dict([(n,k) for n,k in enumerate(clusters.keys())])

matACC_tomatID = dict([(row.rstrip().split('\t')[1],str(row.rstrip().split('\t')[0].split(',')[0])) for row in open(mirnafile).readlines()[1:]])

#############


def ids(self=None, IDs = IDs, LABELS=LABELS, fromTaxo=fromTaxo):
    "ids() functions takes an organism identifier (for human can be 9606 or hsa) and returns all the miRNAs related to that organism. If there is no argument, it returns all the miRNAs in the library."
    
    if not(self):
        to_return = IDs
    else:
        
        try:
            int(self)
            
            if self in list([int(t) for t in obiTaxonomy.common_taxids()]):
                
                to_return = [e for e in IDs if e.split('-')[0]==fromTaxo[self]]                
                
            else:
                to_return = str(self)+': Integer not found'
                
        except ValueError:
            if self in fromTaxo.values():
                
                to_return = [e for e in IDs if e.split('-')[0]==self]
                
            else:
                to_return = self+': String not found'
            
        
    return to_return



class mat_miRNA:
    pass

class pre_miRNA:
    pass
 
      
def miRNA(self, IDs = IDs):
        "miRNA() function takes a miRNA identifier as input and returns a miRNA object."
        if self in IDs:
            attr = [line.rstrip() for line in open('miRNA.txt').readlines()][0].split('\t')
            
            to_return = mat_miRNA()
            setattr(to_return, attr[0], self)
            
            for n,a in enumerate(attr[1:]):
                setattr(to_return, a, miRNA_lib[self][n])            
            
        else:
            to_return = self + ' not found in the miRNA_id library'
        
        return to_return 
    

def premiRNA(self, preIDs = preIDs):
        "premiRNA() function takes a premiRNA identifier as input and returns a premiRNA object."
        if self in preIDs:
            attr = [line.rstrip() for line in open('premiRNA.txt').readlines()][0].split('\t')
            
            to_return = pre_miRNA()
            setattr(to_return, attr[0], self)
            
            for n,a in enumerate(attr[1:]):
                setattr(to_return, a, premiRNA_lib[self][n])            
            
        else:
            to_return = self + ' not found in the premiRNA_id library'
        
        return to_return
    


def get_miRNA_clusters(self, type = 'pre', mat_toPre = mat_toPre, preIDs = preIDs, clusters = clusters, num_clusters=num_clusters):
    "get_miRNA_cluster() function take a mature miRNA or a pre-miRNA identifier and return an integer cluster identifier."
    if type == 'mat':
        print self, 'found in IDs'        
        
        to_return = [(pre1, clusters[pre1]) for pre1 in mat_toPre[self]]
                
    elif type == 'pre':
        print self, 'found in preIDs', 
        to_return = clusters[self]
        
    else:
        to_return = self + ' not found in cluster library'
        
    return to_return


def cluster(self, clusters = clusters, num_clusters = num_clusters):
    "cluster() function take a cluster identifier and return the relative list of clusters."
    try:
        int(self)
        
        to_return = clusters[num_clusters[self]]
    except ValueError:
        
        to_return = self + ' must be an integer for indexing a cluster.'
        
    return to_return


def get_matID(self, matACC_tomatID=matACC_tomatID):
    
    accs = self.split(',')
    
    to_return = []
    for s in accs:
        if s in matACC_tomatID:
            to_return.append(matACC_tomatID[s])
        else:
            to_return.append('None')
            
    return to_return
    
           