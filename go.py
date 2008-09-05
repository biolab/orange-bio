"""A library for handling gene ontologies.
"""
import _GOLib
import re
from sets import Set
from urllib import urlretrieve
from collections import defaultdict
import cPickle
import os

try:
    import orngOrangeFoldersQt4
    default_database_path = orngOrangeFoldersQt4.__getDirectoryNames()['bufferDir']
except:
    default_database_path = os.path.join((os.path.split(__file__)[0] or "."), "data")

data_dir = default_database_path
   
#Currently loaded ontology (GeneOntologyDB)
loadedGO=None

#Currently loaded slim ontology (GeneOntologyDB)
loadedSlimGO=None

#Currently loaded annotation (AnnotationDB)
loadedAnnotation=None

#A dictionary for mapping gene aliases
geneMapper={}

#A dictionary for mapping GOIs's
termMapper={}

namespaceDict={
    "biological_process":1, "P":1,
    "cellular_component":2, "C":2,
    "molecular_function":4, "F":4}

evidenceTypes = {
##Experimental
    'EXP': 'Inferred from Experiment',
    'IDA': 'Inferred from Direct Assay',
    'IPI': 'Inferred from Physical Interaction', ## [with <database:protein_name>]',
    'IMP': 'Inferred from Mutant Phenotype',
    'IGI': 'Inferred from Genetic Interaction', ## [with <database:gene_symbol[allele_symbol]>]',
    'IEP': 'Inferred from Expression Pattern',
##Computational Analysis Evidence Codes
    'ISS': 'Inferred from Sequence Similarity', ## [with <database:sequence_id>] ',
    'ISA': 'Inferred from Sequence Alignment',
    'ISO': 'Inferred from Sequence Orthology',
    'ISM': 'Inferred from Sequence Model',
    'IGC': 'Inferred from Genomic Context',
    'RCA': 'Inferred from Reviewed Computational Analysis',
##Author Statement Evidence Codes
    'TAS': 'Traceable author statement',
    'NAS': 'Non-traceable author statement',
##Curatorial Statement Evidence Codes
    'IC': 'Inferred by curator',
    'ND': 'No biological data available',
##Computationally-assigned Evidence Codes
    'IEA': 'Inferred from electronic annotation', ## [to <database:id>]',
##Obsolete Evidence Codes
    'NR': 'Not Recorded(Obsolete)'
}
##evidenceDict={"IMP":1, "IGI":2, "IPI":4, "ISS":8, "IDA":16, "IEP":32, "IEA":64,
##              "TAS":128, "NAS":256, "ND":512, "IC":1024, "RCA":2048, "IGC":4096, "RCA":8192, "NR":16384}

evidenceDict=defaultdict(int, [(e, 2**i) for i, e in enumerate(evidenceTypes.keys())])

evidenceTypesOrdered = [
'EXP',
'IDA',
'IPI',
'IMP',
'IGI',
'IEP',
##Computational Analysis Evidence Codes
'ISS',
'ISA',
'ISO',
'ISM',
'IGC',
'RCA',
##Author Statement Evidence Codes
'TAS',
'NAS',
##Curatorial Statement Evidence Codes
'IC',
'ND',
##Computationally-assigned Evidence Codes
'IEA',
##Obsolete Evidence Codes
'NR'
]

multiplicitySet=Set(["alt_id","is_a","subset","synonym","related_synonym","exact_synonym","broad_synonym","narrow_synonym",
                     "xref_analog","xref_unknown","relationship"])

annotationFields=["DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GOID", "DB_Reference","Evidence","With_From","Aspect",
                  "DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxon","Date","Assigned_by"]

annotationFieldsDict={"DB":0,
                      "DB_Object_ID":1,
                      "DB_Object_Symbol":2,
                      "Qualifier":3,
                      "GO_ID":4,
                      "GO ID":4,
                      "DB_Reference":5,
                      "DB:Reference":5,
                      "Evidence_code":6,
                      "Evidence code":6,
                      "With_or_From":7,
                      "With (or) From":7,
                      "Aspect":8,
                      "DB_Object_Name":9,
                      "DB_Object_Synonym":10,
                      "DB_Object_Type":11,
                      "taxon":12,
                      "Date":13,
                      "Assigned_by":14}

def __evidenceToInt(evidenceList):
    if not evidenceList:
        return 4095
    evidence=0
    for e in evidenceList:
        if type(e)==str:
            evidence|=evidenceDict[e]
        elif type(e)==int:
            evidence|=e
    return evidence

def __evidenceToList(evidenceCode):
    evidence=[]
    for key, val in evidenceDict.items():
        if val&evidenceCode:
            evidence.append(key)
    return evidence

class AnnotationDB(object):
    """An object holding the annotation database.
    members:
        -geneNames      -- Names of all the genes in the annotation
        -annotationList -- List of instances of Annotation class holding the details for each annotation record
        -aliasMapper    -- Alias mapper maps known aliases to gene names (column3 DB_Object_Symbol of annotation file)
        -geneNamesDict  -- A dictionary mapping any known alias to a list [DB_Object_ID, DB_Object_Symbol [,DB_Object_Synonym]] i.d. all known names
        -geneAnnotations-- A dictionary mapping gene name (DB_Object_Symbol) to a list of all instances of Annotation with this name
    """
    __annotation=None       #holds a C annotation structure for internal use(do not touch!)
    geneNames=None          #a list of all gene names in annotation
    annotationList=None     #a list of all instances of class Annotation
    aliasMapper=None        #alias mapper maps known aliases to gene names (column3 DB_Object_Symbol of annotation file)
    geneNamesDict=None      #a dictionary mapping any known alias to a list [DB_Object_ID, DB_Object_Symbol [,DB_Object_Synonym]] i.e. all known names
    geneAnnotations=None    #a dictionary mapping gene name (DB_Object_Symbol) to a list of all instances of Annotation with this name

class GeneOntologyDB(object):
    """An object holding the ontology database.
    members:
        -terms              - List of instances of GOTerm class holding the details for each term
        -termDict           - A dictionary mapping GO id's to an instance of GOTerm class with that id
        -termDescriptorDict - A dictionary mapping GO id's and alt_id's to a tuple (id, namespace, def, alt_id)
        -aliasMapper        - A dictionary mapping alt_id's and id's to id's
    """
    __ontology=None         #holds a C ontology structure for internal use(do not touch!)
    terms=None              #holds a list of all instances of class GOTerm
    termDict=None           #a dictionary mapping GO id's to instances of class GOTerm
    termDescriptorDict=None #a dictionary mapping GO id's and alt_id's to a tuple (id, namespace, def, alt_id)
    aliasMapper=None        #a dictionary mapping alt_id's and id's to id's

class GOTerm(object):
    """Holds the data for one term read from the ontology file. All fields from the ontology can be accsessed by their
    original name (e.g. the is_a field in the ontology can be accsessed like term.is_a), except for the def field that
    interferes with the python def statment and can be accesed like term.def_ or term.__dict__['def'].
    The fields that can have multiple values are stored as lists of strings otherwise the string after the field name from the
    ontology is supplied. If a field is not defined accessing it will result in an exception.
    The object also provides the folowing data memebers for quicker accsess: GOId, aspect, parents(list of GOId's of parents terms)."""
    def __init__(self, term):
        self.original={}
        self.GOId=None
        self.parents=[]
        self.rType={} #relationship type with a parent
        self.fullText=""
        self.aspect="unknown"
        self.alt_id=[]

        self.processBlock(term)

    def processBlock(self, text):
        for s in text.split("\n"):
            if ":" not in s:
                continue
            index=s.find(":")
            block=s[:index]
            body=s[index+1 :].strip(" ")
            if block in multiplicitySet:
                if self.__dict__.has_key(block):
                    self.__dict__[block].append(body)
                else:
                    self.__dict__[block]=[body]
            else:
                self.__dict__[block]=body            
        self.GOId=self.id
        self.parents=self.__dict__.get("is_a", [])
        self.parents=[s.split("!")[0].strip(" ") for s in self.parents]
        for p in self.parents:
            self.rType[p]="is_a"
        for r in self.__dict__.get("relationship",[]):
            s=r.split("!")[0].strip(" ").split(" ")
            rT,id=s[0],s[1]
            self.parents.append(id)
            self.rType[id]=rT
        self.aspect=self.__dict__.get("namespace", "unknown")
        self.alt_id=self.__dict__.get("alt_id", [])
        self.fullText+=text

    def __getattr__(self, name):
        if name=="def_":
            return self.__dict__["def"]
        else:
            raise AttributeError(name)
        
    def toTuple(self):
        return (self.GOId, self.parents)
    
_re_obj_name_ = re.compile("([a-zA-z0-9-_]+)")

class Annotation(object):
    """Holds the data for an annotation record read from the annotation file. Fields can be
    accessed with the names: DB, DB_Object_ID, DB_Object_Symbol, Qualifier, GO_ID, DB_Reference,
    Evidence_code, With_or_From, Aspect, DB_Object_Name, DB_Object_Synonym, DB_Object_Type, taxon,
    Date, Assigned_by (e.g. rec.GO_ID)
    or by supplying the original name of the field (see http://geneontology.org/GO.annotation.shtml#file)
    to the get method (e.g. rec.get("GO ID"))
    The object also provides the folowing data members for quicker access: geneName, GOId, evidence,
    aspect and alias(a list of aliases)
    """
    def __init__(self, fullText):
        self.fullText=fullText
        self.original=fullText.split("\t")
        self.geneName=self.original[2]
        self.GOId=self.original[4].strip(" ")
        self.evidence=self.original[6].strip(" ")
        self.aspect=self.original[8].strip(" ")
        self.alias=self.original[10].split("|") 
        for key, val in zip(annotationFields, self.original):
            self.__dict__[key]=val

        self.aditionalAliases = []
        if ":" in self.DB_Object_Name:
            self.aditionalAliases = _re_obj_name_.findall(self.DB_Object_Name.split(":")[0])

    def __getattr__(self, name):
        if name in annotationFieldsDict:
            return self.original[annotationFieldsDict[name]]
        else:
            raise AttributeError(name)

    def get(self, name):
        if name in annotationFieldsDict:
            return self.original[annotationFieldsDict[name]]
        else:
            raise ValueError(name)
        
    def toTuple(self):
        return (self.geneName, self.GOId, evidenceDict[self.evidence] , namespaceDict[self.aspect]) #evidenceMapper[self.evidence]

def setDataDir(datadir):
    """Set the data directory where the annotation and ontology files are stored (by default a directory named data
    located where the GOLib is instaled e.g. ...site-packages/GOLib/data)"""
    global data_dir
    data_dir=datadir

def getDataDir():
    """Get the data directory where the annotation and ontology files are stored (by default a directory named data
    located where the GOLib is instaled e.g. ...site-packages/GOLib/data)"""
    return data_dir

def loadAnnotation(organism="sgd", forceReload=False, progressCallback=None):
    """Loads the annotation for the specified organism"""
    global loadedAnnotation
    if loadedAnnotation and loadedAnnotation.__file__.endswith(organism) and not forceReload:
        return
    loadedAnnotation=loadAnnotationFrom(os.path.join(data_dir,"gene_association."+organism), progressCallback)#+".PyAnnotationDB")
    global geneMapper
    geneMapper=loadedAnnotation.aliasMapper

def loadGO(forceReload=False, progressCallback=None):
    """Loads the ontology from 'data_dir//gene_ontlogy.obo' where data_dir is set using setDataDir (default: .//data)"""
    global loadedGO
    if loadedGO and not forceReload:
        return
    loadedGO=loadOntologyFrom(os.path.join(data_dir,"gene_ontology.obo"), progressCallback)#.PyOntologyDB")
    global termMapper
    termMapper=loadedGO.aliasMapper
    
def mapTermId(TermId):
    """Maps the TermId to id if TermId a known alias for id (TermId can map to it self), if TermId unknown return None""" 
    return loadedGO.aliasMapper.get(TermId, None)

def mapGeneName(GeneName):
    """Maps the GeneName to name if GeneName a known alias for name (GeneName can map to it self), if GeneName unknown return None""" 
    return loadedAnnotation.aliasMapper.get(GeneName, None)

def mapGeneNames(names=[]):
    return filter(lambda a:a, map(mapGeneName, names))

def __filterSlimsGOId(d):
    slims=[t.GOId for t in getSlims()]
    slims=filter(lambda id:id in slims, d.keys())
    return dict([(id, d[id]) for id in slims])

def GOTermFinder(clusterGeneList, referenceGenes=None, evidenceCodes=None, slimsOnly=False, aspect="P", progressCallback=None):
    """The method accepts a list of cluster genes, optionally a list of reference genes (otherwise all annotated genes appearing in the loaded annotation file are used),
    and optionally a list of annotation evidence codes to use, otherwise all evidence codes are used. The slimsOnly argument indicates if only GO slims are to be used,
    otherwise all GO terms are considered. The method returns a dictionary of significant GO terms, items are (cluster genes that map to each selected term, p-value,
    number of reference genes that map to each significant GO term). The progressCallback if present will be called with a single argument for all values in range [0..99]
    """
    if not referenceGenes:
        referenceGenes=loadedAnnotation.geneNames
    if slimsOnly and not loadedSlimGO:
        setSlims()
    #clusterGeneList=mapGeneNames(clusterGeneList)
    #referenceGenes=mapGeneNames(referenceGenes)
        
    annotation=loadedAnnotation.__annotation
    goslim=loadedSlimGO and loadedSlimGO.__ontology
    evidence=__evidenceToInt(evidenceCodes)
    aspect=namespaceDict[aspect]
    result=_GOLib.GOTermFinder(clusterGeneList, referenceGenes, slimsOnly, evidence, aspect, annotation, loadedGO.__ontology, goslim, progressCallback)
    if slimsOnly:
        return __filterSlimsGOId(result)
    else:
        return result
    
def findTerms(geneList, slimsOnly=False, aspect=["F","C","P"], directAnnotationOnly=False, evidenceCodes=None, reportEvidence=True, progressCallback=None):
    """For each gene in geneList search for matching GO terms. Argument slimsOnly restricts the GO terms to the slim set. The method returns a dictionary where key is a
    matching GO term and items are (gene, evidence) if reportEvidence == True [gene only, if reportEvidence=False] that map to the term. Climb the GO if directAnnotationOnly=False,
    otherwise report direct annotation only.
    """
    evidence=__evidenceToInt(evidenceCodes)
    if slimsOnly and not loadedSlimGO:
        setSlims()
    aa=0
    for a in aspect:
        aa|=namespaceDict[a]
    #geneList=mapGeneNames(geneList)
    goslim=loadedSlimGO and loadedSlimGO.__ontology
    annotation=loadedAnnotation.__annotation
    result=_GOLib.findTerms(geneList, slimsOnly, directAnnotationOnly, aa, evidence, reportEvidence, annotation, loadedGO.__ontology, goslim, progressCallback)
    if slimsOnly:
        result=__filterSlimsGOId(result)
    if reportEvidence:
        result=dict([(key, [(gene, __evidenceToList(evidence)) for gene ,evidence in val]) for key, val in result.items()])
    return result

def findGenes(GOTerms=[], directAnnotationOnly=False, evidenceCodes=None, reportEvidence=True, progressCallback=None):
    """Return a dictionary where key is a matching gene and items are (GO terms) or (GO term, list of evidences) from the GOterms list.
    (Note this will take a lot of time if the directAnnotationOnly=False)"""
    evidence=__evidenceToInt(evidenceCodes)
    result=_GOLib.findGenes(GOTerms, evidence, reportEvidence, directAnnotationOnly, loadedAnnotation.__annotation, loadedGO.__ontology, progressCallback)
    if reportEvidence:
        result=dict([(key,[(term, __evidenceToList(evidence)) for term, evidence in val]) for key, val in result.items()])
    return result        
    
def extractGODAG(GOTerms=[]):
    """Return the part of GO DAG that includes the listed GOTerms."""
    expanded=[]
    queue=list(GOTerms)
    while queue:
        id=queue.pop(0)
        term=loadedGO.termDict.get(id, None)
        term=term or loadedGO.termDict.get(loadedGO.aliasMapper.get(id, None), None)
        if term and (term.id not in expanded):
            expanded.append(term.id)
            queue.extend(term.parents)
    terms=[loadedGO.termDict[id] for id in expanded if id in loadedGO.termDict]
    return terms      

def __DAGDepth(term, cache={}):
	if term.parents:
		d=[]
		for t in term.parents:
			if t in cache:
				d.append(cache[t]+1)
			else:
				d.append(__DAGDepth(loadedGO.termDict[t], cache)+1)
		depth=min(d)
		cache[term.id]=depth
		return depth
	else:
		return 1

def DAGDepth(DAGTerms=[]):
    """returns the maximum depth of terms in DAGTerms"""
    cache={}
    #DAGTerms=[type(t)==GOTerm and t or loadedGO.termDict[t] for t in DAGTerms]
    return max([__DAGDepth(t, cache) for t in DAGTerms])

def DAGFilterForDepth(DAGTerms, depth):
    cache={}
    return filter(lambda t:__DAGDepth(t,cache)<=depth, DAGTerms)

def extractGODAGToDisplay(geneList):
    """returns the part of the GO that all of the genes in geneList map to
    """
    terms=findTerms(geneList)
    return terms.keys()
    
def setSlims(slims=None):
    """Set the slim subset of a loaded ontology. slims can be:
        -a string identifier of a subsetdef: (e.g. "goslim_generic", "goslim_plant" ...)
        -a filename of a slim ontology (e.g. "goslim_generic.obo" ...)
        -a list of GO term id's
    """
    global loadedSlimGO
    goslims=[]
    slims= slims or "goslim_generic"
    try:
        loadedSlimGO=loadOntologyFrom(slims)
    except:
        if type(slims)==list:
            goslims=[loadedGO.termDict[term] for term in slims]
        else:
            goslims=filter(lambda t:slims in t.__dict__.get("subset",[]) , loadedGO.terms)
        loadedSlimGO=GeneOntologyDB()
        loadedSlimGO.__ontology=_GOLib.parseGOTerms([g.toTuple() for g in goslims], loadedGO.aliasMapper)
        #loadedSlimGO.__ontology.aliasMapper=loadedGO.aliasMapper
        loadedSlimGO.terms=goslims
        loadedSlimGO.termDict=loadedGO.termDict
        loadedSlimGO.aliasMapper=loadedGO.aliasMapper
        loadedSlimGO.termDescriptorDict=loadedGO.termDescriptorDict

def getSlims():
    """Returns the curently loaded slim terms"""
    return loadedSlimGO.terms

class __progressCallWrapper:
        def __init__(self,callback):
            self.callback=callback
        def __call__(self, bCount, bSize, fSize):
            #print bCount, bSize, fSize
            if fSize==-1:
                fSize=10000000
            self.callback(100*bCount*bSize/fSize)
            
def downloadGO(progressCallback=None):
    """Downloads the curent gene ontology from http://www.geneontology.org/ontology/gene_ontology.obo"""
    downloadGOTo(progressCallback=None)
##    urlretrieve("http://www.geneontology.org/ontology/gene_ontology.obo", os.path.join(data_dir, "gene_ontology.obo"), progressCallback and __progressCallWrapper(progressCallback))
##    file=open(os.path.join(data_dir,"gene_ontology.obo"))
##    data=file.read()
##    c=re.compile("\[Term\].*?\n\n",re.DOTALL)
##    match=c.findall(data)
##    go=parseGeneOntology(match)
    #cPickle.dump(go, open(data_dir+"gene_ontology.obo.PyOntologyDB", "w"))

def downloadGOTo(filename=None, progressCallback=None):
    """Downloads the curent gene ontology from http://www.geneontology.org/ontology/gene_ontology.obo, to filename"""
    filename=filename or os.path.join(data_dir,"gene_ontology.obo")
    urlretrieve("http://www.geneontology.org/ontology/gene_ontology.obo", filename, progressCallback and __progressCallWrapper(progressCallback))
    file=open(filename)
    data=file.read()
    c=re.compile("\[Term\].*?\n\n",re.DOTALL)
    match=c.findall(data)
    go=parseGeneOntology(match)
    #cPickle.dump(go, open(filename+".PyOntologyDB", "w"))
    
def downloadAnnotation(organism="sgd", progressCallback=None):
    """Downloads the annotation for the specified organism (e.g. "sgd", "fb", "mgi",...)"""
        
    #urlretrieve("http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association."+organism+".gz",
    #            data_dir+"//gene_association."+organism+".gz", progressCallback and __progressCallWrapper(progressCallback))
    urlretrieve("http://www.geneontology.org/gene-associations/gene_association."+organism+".gz",
                os.path.join(data_dir,"gene_association."+organism+".gz"), progressCallback and __progressCallWrapper(progressCallback))
    from gzip import GzipFile
    gfile=GzipFile(os.path.join(data_dir, "gene_association."+organism+".gz"),"r")
    data=gfile.readlines()
    file=open(os.path.join(data_dir, "gene_association."+organism),"w")
    file.writelines(data)
    #__splitAnnotation(data, organism)
    anno=parseAnnotation(data)
    import cPickle
    cPickle.dump(anno.aliasMapper.keys(), open(os.path.join(data_dir, "gene_names."+organism), "w"))

def downloadAnnotationTo(organism="sgd", filename=None, progressCallback=None):
    filename = filename or os.path.join(data_dir, "gene_association."+organism+".gz")
    urlretrieve("http://www.geneontology.org/gene-associations/gene_association."+organism+".gz",
                filename+".gz", progressCallback and __progressCallWrapper(progressCallback))
    from gzip import GzipFile
    gfile=GzipFile(filename+".gz", "r")
    data=gfile.readlines()
    file=open(filename,"w")
    file.writelines(data)
    #__splitAnnotation(data, organism)
    anno=parseAnnotation(data)
    import cPickle
    cPickle.dump(anno.aliasMapper.keys(), open(os.path.join((os.path.split(filename)[0] or "."), "gene_names."+organism), "w"))

def getCachedGeneNames(organism="sgd"):
    import cPickle
    return cPickle.load(open(os.path.join(data_dir, "gene_names."+organism)))

def listOrganisms():
    """Connect to http://www.geneontology.org/GO.current.annotations.shtml, parse out the organism names
    appearing in the table, and return the list of organisms."""
    try:
        urlretrieve("http://www.geneontology.org/GO.current.annotations.shtml", os.path.join(data_dir, "annotations.shtml"))
    except:
        print "Failed to connect to http://www.geneontology.org/GO.current.annotations.shtml. Trying to find a local copy"
    file=open(os.path.join(data_dir, "annotations.shtml"))
    data=file.read()
    #match=re.findall(r'http://www\.geneontology\.org/cgi-bin/downloadGOGA\.pl/gene_association\..+?gz', data)
    match=re.findall(r'http://cvsweb\.geneontology\.org/cgi-bin/cvsweb\.cgi/go/gene-associations/gene_association\..+?gz\?rev=HEAD', data)
    organisms=[]
    for m in match:
        #print m
        #names=re.findall(r'>.+?<', m)
        organisms.append(m.split(".")[-2])
    return organisms

def listDownloadedOrganisms():
    """Returns a list with organism names off all local downloaded annotations"""
    import os
    files=os.listdir(data_dir)
    files=filter(lambda n: n.startswith("gene_association") and not n.endswith(".gz"), files)
    return [s[17:] for s in files]

def parseGeneOntology(data, progressCallback=None):
    terms=[]
    termDict={}
    aliasMapper={}
    goTermDict={}
    datalen=len(data)
    milestones=Set(range(0,datalen,max(datalen/100,1)))
    for i, term in enumerate(data):
        t=GOTerm(term)
        termDict[t.id]=t
        terms.append(t)
        alt=t.__dict__.get("alt_id", [])
        aliasMapper.update(dict([(alt_id.strip(" "), t.id) for alt_id in alt]))
        aliasMapper[t.id]=t.id
        tt=(t.id, t.namespace, t.__dict__.get("def",""), alt)
        for alt_id in alt:
            goTermDict[alt_id]=tt
        goTermDict[t.id]=tt
        if progressCallback and i in milestones:
            progressCallback(100.0*i/datalen)
    GO=GeneOntologyDB()
    GO.terms=terms
    GO.termDict=termDict
    GO.aliasMapper=aliasMapper
    GO.termDescriptorDict=goTermDict
    return GO

def loadOntologyFrom(filename,progressCallback=None):
    if filename.endswith(".PyOntologyDB"):
        db=cPickle.load(open(filename))
    else:
        file=open(filename)
        data=file.read()
        c=re.compile("\[Term\].*?\n\n",re.DOTALL)
        match=c.findall(data)
        db=parseGeneOntology(match, progressCallback)
    db.__ontology=_GOLib.parseGOTerms([t.toTuple() for t in db.terms], db.aliasMapper)
    #db.__ontology.aliasMapper=db.aliasMapper
    db.__file__=filename
    return db

def parseAnnotation(data, progressCallback=None):
    aliasMapper={}
    annotationList=[]
    geneNames=Set()
    geneNamesDict={}
    geneAnnotation={}
    datalen=len(data)
    milestones=Set(range(0,datalen, max(datalen/100,1)))
    for i,line in enumerate(data):
        if line.startswith("!"):
            continue
        a=Annotation(line)
        if not a.geneName or not a.GOId:
            continue
        if a.geneName not in geneNames:
            geneNames.add(a.geneName)
            geneAnnotation[a.geneName]=[a]
            for alias in a.alias:
                aliasMapper[alias]=a.geneName
            for alias in a.aditionalAliases:
                aliasMapper[alias]=a.geneName
            aliasMapper[a.geneName]=a.geneName
            aliasMapper[a.DB_Object_ID]=a.geneName
            names=[a.original[1], a.original[2]]
            names.extend(a.alias)
            for n in names:
                geneNamesDict[n]=names
        else:
            geneAnnotation[a.geneName].append(a)
        annotationList.append(a)
        if progressCallback and i in milestones:
            progressCallback(100.0*i/datalen)
    a=AnnotationDB()
    a.annotationList=annotationList
    a.aliasMapper=aliasMapper
    a.geneNames=list(geneNames)
    a.geneNamesDict=geneNamesDict
    a.geneAnnotations=geneAnnotation
    return a

def loadAnnotationFrom(filename, progressCallback=None):
    if filename.endswith(".PyAnnotationDB"):
        anno=cPickle.load(open(filename))
    else:
        file=open(filename)
        data=file.readlines()
        anno=parseAnnotation(data, progressCallback)

    aa=[a.toTuple() for a in anno.annotationList if "NOT" not in a.Qualifier]
    aa.sort()
    anno.__annotation=_GOLib.parseAnnotation(aa, anno.aliasMapper)
    #anno.__annotation.aliasMapper=anno.aliasMapper
    anno.__file__=filename
    return anno

def filterByPValue(terms, maxPValue=0.1):
    """Filters the terms by the p-value. Asumes terms is is a dict with the same structure as returned from GOTermFinderFunc
    """
    return dict(filter(lambda (k,e): e[1]<maxPValue, terms.items()))

def filterByFrequency(terms, minF=2):
    """Filters the terms by the cluster frequency. Asumes terms is is a dict with the same structure as returned from GOTermFinderFunc
    """
    return dict(filter(lambda (k,e): len(e[0])>=minF, terms.items()))

def filterByRefFrequency(terms, minF=4):
    """Filters the terms by the reference frequency. Asumes terms is is a dict with the same structure as returned from GOTermFinderFunc
    """
    return dict(filter(lambda (k,e): e[2]>=minF, terms.items()))

def drawEnrichmentGraph(termsList, clusterSize, refSize, filename="graph.png", width=None, height=None):
	drawEnrichmentGraph_tostream(termsList, clusterSize, refSize, open(filename, "wb"), width, height)

def drawEnrichmentGraph_tostream(GOTerms, clusterSize, refSize, fh, width=None, height=None):
    def getParents(term):
        parents = extractGODAG([term])
        parents = filter(lambda t: t.id in GOTerms and t.id!=term, parents)
        c = []
        map(c.extend, [getParents(t.id) for t in parents])
        parents = filter(lambda t: t not in c, parents)
        return parents
    parents = dict([(term, getParents(term)) for term in GOTerms])
    #print "Parentes", parents
    def getChildren(term):
        return filter(lambda t: term in [p.id for p in parents[t]], GOTerms.keys())
    topLevelTerms = filter(lambda t: not parents[t], parents.keys())
    #print "Top level terms", topLevelTerms
    termsList=[]
    def collect(term, parent):
        termsList.append(
            ((float(len(GOTerms[term][0]))/clusterSize) / (float(GOTerms[term][2])/refSize),
            len(GOTerms[term][0]),
            GOTerms[term][2],
            "%.4f" % GOTerms[term][1],
            loadedGO.termDict[term].name,
            loadedGO.termDict[term].id,
            ", ".join(GOTerms[term][0]),
            parent)
            )
##        print float(len(GOTerms[term][0])), float(GOTerms[term][2]), clusterSize, refSize
        parent = len(termsList)-1
        for c in getChildren(term):
            collect(c, parent)
                         
    for topTerm in topLevelTerms:
        collect(topTerm, None)

    drawEnrichmentGraphPIL_tostream(termsList, fh, width, height)
        
def drawEnrichmentGraphPIL_tostream(termsList, fh, width=None, height=None):
    from PIL import Image, ImageDraw, ImageFont
    backgroundColor = (255, 255, 255)
    textColor = (0, 0, 0)
    graphColor = (0, 0, 255)
    fontSize = height==None and 12 or (height-60)/len(termsList)
    font = ImageFont.load_default()
    try:
        font = ImageFont.truetype("arial.ttf", fontSize)
    except:
        pass
    getMaxTextHeightHint = lambda l: max([font.getsize(t)[1] for t in l])
    getMaxTextWidthHint = lambda l: max([font.getsize(t)[0] for t in l])
    maxFoldWidth = width!=None and min(150, width/6) or 150
    maxFoldEnrichment = max([t[0] for t in termsList])
    foldNormalizationFactor = float(maxFoldWidth)/maxFoldEnrichment
    foldWidths = [int(foldNormalizationFactor*term[0]) for term in termsList]
    treeStep = 10
    treeWidth = {}
    for i, term in enumerate(termsList):
        treeWidth[i] = (term[7]==None and 1 or treeWidth[term[7]]+1)
    treeStep = width!=None and min(treeStep, width/(6*max(treeWidth.values())) or 2) or treeStep
    treeWidth = [w*treeStep + foldWidths[i] for i, w in treeWidth.items()]
    treeWidth = max(treeWidth) - maxFoldWidth
    verticalMargin = 10
    horizontalMargin = 10
##    print verticalMargin, maxFoldWidth, treeWidth
##    treeWidth = 100
    firstColumnStart = verticalMargin + maxFoldWidth + treeWidth + 10
    secondColumnStart = firstColumnStart + getMaxTextWidthHint([str(t[1]) for t in termsList]+["List"]) + 2
    thirdColumnStart = secondColumnStart + getMaxTextWidthHint([str(t[2]) for t in termsList]+["Total"]) + 2
    fourthColumnStart = thirdColumnStart + getMaxTextWidthHint([str(t[3]) for t in termsList]+["p-value"]) + 4
##    maxAnnotationTextWidth = width==None and getMaxTextWidthHint([str(t[4]) for t in termsList]+["Annotation"]) or (width - fourthColumnStart - verticalMargin) * 2 / 3
    maxAnnotationTextWidth = width==None and getMaxTextWidthHint([str(t[4]) for t in termsList]+["Annotation"]) or max((width - fourthColumnStart - verticalMargin) * 2 / 3, getMaxTextWidthHint([t[5] for t in termsList]+["Annotation"]))
    fifthColumnStart  = fourthColumnStart + maxAnnotationTextWidth + 4
    maxGenesTextWidth = width==None and getMaxTextWidthHint([str(t[5]) for t in termsList]+["Genes"]) or (width - fourthColumnStart - verticalMargin) / 3
    
    legendHeight = font.getsize("1234567890")[1]*2
    termHeight = font.getsize("A")[1]
##    print fourthColumnStart, maxAnnotationTextWidth, verticalMargin
    width = fifthColumnStart + maxGenesTextWidth + verticalMargin
    height = len(termsList)*termHeight+2*(legendHeight+horizontalMargin)

    image = Image.new("RGB", (width, height), backgroundColor)
    draw = ImageDraw.Draw(image)

    def truncText(text, maxWidth, append=""):
        #print getMaxTextWidthHint([text]), maxAnnotationTextWidth
        if getMaxTextWidthHint([text])>maxWidth:
            while getMaxTextWidthHint([text+"..."+append])>maxWidth and text:
                text = text[:-1]
            if text:
                text = text+"..."+append
            else:
                text = append
        return text
    currentY = horizontalMargin + legendHeight
    connectAtX = {}
    for i, term in enumerate(termsList):
        draw.line([(verticalMargin, currentY+termHeight/2), (verticalMargin + foldWidths[i], currentY+termHeight/2)], width=termHeight-2, fill=graphColor)
        draw.text((firstColumnStart, currentY), str(term[1]), font=font, fill=textColor)
        draw.text((secondColumnStart, currentY), str(term[2]), font=font, fill=textColor)
        draw.text((thirdColumnStart, currentY), str(term[3]), font=font, fill=textColor)
        annotText = width!=None and truncText(str(term[4]), maxAnnotationTextWidth, str(term[5])) or str(term[4])
        draw.text((fourthColumnStart, currentY), annotText, font=font, fill=textColor)
        genesText = width!=None and truncText(str(term[6]), maxGenesTextWidth) or str(term[6])
        draw.text((fifthColumnStart, currentY), genesText, font=font, fill=textColor)
        lineEnd = term[7]==None and firstColumnStart-10 or connectAtX[term[7]]
        draw.line([(verticalMargin+foldWidths[i]+1, currentY+termHeight/2), (lineEnd, currentY+termHeight/2)], width=1, fill=textColor)
        if term[7]!=None:
            draw.line([(lineEnd, currentY+termHeight/2), (lineEnd, currentY+termHeight/2 - termHeight*(i-term[7]))], width=1, fill=textColor)
        connectAtX[i] = lineEnd - treeStep
        currentY+=termHeight

    currentY = horizontalMargin
    draw.text((firstColumnStart, currentY), "list", font=font, fill=textColor)
    draw.text((secondColumnStart, currentY), "total", font=font, fill=textColor)
    draw.text((thirdColumnStart, currentY), "p-value", font=font, fill=textColor)
    draw.text((fourthColumnStart, currentY), "Annnotation", font=font, fill=textColor)
    draw.text((fifthColumnStart, currentY), "Genes", font=font, fill=textColor)

    horizontalMargin = 0
    #draw.line([(verticalMargin, height - horizontalMargin - legendHeight), (verticalMargin + maxFoldWidth, height - horizontalMargin - legendHeight)], width=1, fill=textColor)
    draw.line([(verticalMargin, horizontalMargin + legendHeight), (verticalMargin + maxFoldWidth, horizontalMargin + legendHeight)], width=1, fill=textColor)
    maxLabelWidth = getMaxTextWidthHint([" "+str(i) for i in range(int(maxFoldEnrichment+1))])
    numOfLegendLabels = max(int(maxFoldWidth/maxLabelWidth), 2)
    for i in range(numOfLegendLabels+1):
        #draw.line([(verticalMargin + i*maxFoldWidth/10, height - horizontalMargin - legendHeight/2), (verticalMargin + i*maxFoldWidth/10, height - horizontalMargin - legendHeight)], width=1, fill=textColor)
        #draw.text((verticalMargin + i*maxFoldWidth/10 - font.getsize(str(i))[0]/2, height - horizontalMargin - legendHeight/2), str(i), font=font, fill=textColor)

        label = str(int(i*maxFoldEnrichment/numOfLegendLabels))
        draw.line([(verticalMargin + i*maxFoldWidth/numOfLegendLabels, horizontalMargin + legendHeight/2), (verticalMargin + i*maxFoldWidth/numOfLegendLabels, horizontalMargin + legendHeight)], width=1, fill=textColor)
        draw.text((verticalMargin + i*maxFoldWidth/numOfLegendLabels - font.getsize(label)[0]/2, horizontalMargin), label, font=font, fill=textColor)
        
    image.save(fh)

def __test1():
    setDataDir("E://orangecvs//GOLib//data")
    print "Loading GO"
    loadGO()
    print "Loading annotation"
    loadAnnotation()
##    print len(loadedAnnotation.geneNames)
    terms = GOTermFinder(loadedAnnotation.geneNames[:30], aspect="P")
    terms = filterByPValue(terms, 0.05)
    terms = filterByFrequency(terms, 3)
    terms = filterByRefFrequency(terms, 10)
    print terms
    try:
        drawEnrichmentGraph(terms, 30, len(loadedAnnotation.geneNames), filename="pict.png", width=400, height=2000) #, width=400)#, height=1000)
    except Exception, err:
        print err
        raw_input()
    
def __test():
    def call(i): print i
    setDataDir("E://orangecvs//GOLib//data")
    print "Loading GO"
    loadGO()
    print "Loading annotation"
    loadAnnotation()
    print "Loading slims"
    setSlims()
    print "Finding terms"
    print GOTermFinder(loadedAnnotation.geneNames[:20], progressCallback=call)
    print "Finding slim terms"
    terms = GOTermFinder(loadedAnnotation.geneNames[20:30], slimsOnly=True, aspect="P")
    print terms
    drawGOAT(terms)
    print "Finding terms"
    print findTerms(loadedAnnotation.geneNames[100:120])
    print "Finding direct slim terms"
    terms=findTerms(loadedAnnotation.geneNames[200:210], slimsOnly=True, directAnnotationOnly=True)
    print terms
    print "Extracting GO Dag"
    print extractGODAG(terms.keys())
    print "Finding genes"
    print findGenes(terms.keys()[:min(len(terms.keys()),3)], progressCallback=call)#,directAnnotationOnly=True)
    print findGenes(["GO:0005763","GO:0003735","GO:0042255","GO:0043037"], progressCallback=call)#,directAnnotationOnly=True)

if not listDownloadedOrganisms():
    print "Warning!!! No downloaded annotations found!!!"
    print "You can download annotations using the downloadAnnotation function."
    print "e.g. go.downloadAnnotation('sgd')"
try:
    open(os.path.join(data_dir, "gene_ontology.obo"))
except:
    print "Warning!!! No downloaded ontology found!!!"
    print "You can download it using the downloadGO function."
    print "e.g. go.downloadGO()"
    
    
if __name__=="__main__":
    __test1()
    
    
