"""obiGO is a Gene Ontology (GO) Handling Library.

"""
from go import *

from urllib import urlretrieve
from gzip import GzipFile
import tarfile
import shutil
import sys, os
from datetime import datetime
from collections import defaultdict

try:
    import orngServerFiles
    default_database_path = os.path.join(orngServerFiles.localpath(), "GO")
except Exception:
    default_database_path = os.curdir

builtinOBOObjects = ["""
[Typedef]
id: is_a
name: is_a
range: OBO:TERM_OR_TYPE
domain: OBO:TERM_OR_TYPE
definition: The basic subclassing relationship [OBO:defs]"""
,
"""[Typedef]
id: disjoint_from
name: disjoint_from
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that two classes are disjoint [OBO:defs]"""
,
"""[Typedef]
id: instance_of
name: instance_of
range: OBO:TERM
domain: OBO:INSTANCE
definition: Indicates the type of an instance [OBO:defs]"""
,
"""[Typedef]
id: inverse_of
name: inverse_of
range: OBO:TYPE
domain: OBO:TYPE
definition: Indicates that one relationship type is the inverse of another [OBO:defs]"""
,
"""[Typedef]
id: union_of
name: union_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the union of several others [OBO:defs]"""
,
"""[Typedef]
id: intersection_of
name: intersection_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the intersection of several others [OBO:defs]"""]

class OBOObject(object):
    """ Represents a generic OBO object (e.g. Term, Typedef, Instance, ...)
    Example:
    >>> OBOObject(r"[Term]\nid: FOO:001\nname: bar", ontology)
    """
    def __init__(self, stanza=None, ontology=None):
        self.ontology = ontology
        self._lines = []
        self.values = {}
        self.related = set()
        self.relatedTo = set()
        if stanza:
            self.ParseStanza(stanza)

    def ParseStanza(self, stanza):
        for line in stanza.split("\n"):
            if ":" not in line:
                continue
            tag, rest = line.split(":", 1)
            value, modifiers, comment = "", "", ""
            if "!" in rest:
                rest, comment = rest.split("!")
            if "{" in rest:
                value, modifiers = rest.split("{", 1)
                modifiers = modifiers.strip("}")
            else:
                value = rest
            value = value.strip()
            self._lines.append((tag, value, modifiers, comment))
            if tag in multipleTagSet:
                self.values[tag] = self.values.get(tag, []) + [value]
            else:
                self.values[tag] = value
        self.related = set(self.GetRelatedObjects())

    def GetRelatedObjects(self):
        """ Return a list of tuple pairs where the first element is relationship
        typeId and the second id of object to whom the relationship applys to.
        """
        result = [(typeId, id) for typeId in ["is_a"] for id in self.values.get(typeId, [])] ##TODO add other defined Typedef ids
        result = result + [tuple(r.split(None, 1)) for r in self.values.get("relationship", [])]
        return result

    def __repr__(self):
        """ Return a string representation of the object in OBO format
        """
        repr = "[%s]\n" % type(self).__name__
        for tag, value, modifiers, comment in self._lines:
            repr = repr + tag + ": " + value
            if modifiers:
                repr = repr + "{ " + modifiers + " }"
            if comment:
                repr = repr + " ! " + comment
            repr = repr + "\n"
        return repr

    def __str__(self):
        """ Return the OBO object id entry
        """
        return "%s: %s" % (self.id, self.name)

    def __getattr__(self, tag):
        """ Return value for the tag
        """
        try:
            if tag == "def_":
                return self.values["def"]
            else:
                return self.values[tag]
        except KeyError:
            raise AttributeError(tag)

    def __iter__(self):
        """ Iterates over sub terms
        """
        for typeId, id in self.relatedTo:
            yield (typeId, self.ontology[id])
        
class Term(OBOObject):
    pass

class Typedef(OBOObject):
    pass

class Instance(OBOObject):
    pass
        
class Ontology(object):
    """Ontology is the main class representing a gene ontology."""
    def __init__(self, file=None, progressCallback=None):
        """ Initialize the ontology from file. The optional progressCallback will be called with a single argument to report on the progress.
        """
        self.terms = {}
        self.typedefs = {}
        self.instances = {}
        self.slimsSubset = set()
        if file and os.path.exists(file):
            self.ParseFile(file, progressCallback)
        else:
            fool = self.Load(progressCallback)
            self.__dict__ = fool.__dict__ ## A fool and his attributes are soon parted

    @classmethod
    def Load(cls, progressCallback=None):
        """ A class method that tries to load the ontology file from default_database_path. It looks for a filename starting with 'gene_ontology'.
        """
        filename = os.path.join(default_database_path, "gene_ontology_edit.obo.tar.gz")
        if not os.path.isfile(filename) and not os.path.isdir(filename):
##            print "Ontology file not found on local disk"
##            print "Downloading ontology ..."
            import orngServerFiles
            orngServerFiles.download("GO", "gene_ontology_edit.obo.tar.gz")
        try:
            return cls(filename, progressCallback=progressCallback)
        except (IOError, OSError), ex:
            print ex
            raise Exception, "Could not locate ontology file"
        
    def ParseFile(self, file, progressCallback=None):
        """ Parse the file. file can be a filename string or an open filelike object. The optional progressCallback will be called with a single argument to report on the progress.
        """
        if type(file) == str:
            if os.path.isfile(file) and tarfile.is_tarfile(file):
                f = tarfile.open(file).extractfile("gene_ontology_edit.obo")
            elif os.path.isfile(file):
                f = open(file)
            else:
                f = open(os.path.join(file, "gene_ontology_edit.obo"))
        else:
            f = file
        
        data = f.readlines()
        data = "".join([line for line in data if not line.startswith("!")])
        self.header = data[: data.index("[Term]")]
        c=re.compile("\[.+?\].*?\n\n", re.DOTALL)
        data=c.findall(data)

        milestones = set(i for i in range(0, len(data), max(len(data)/100, 1)))
        for i, block in enumerate(builtinOBOObjects + data):
            if block.startswith("[Term]"):
                term = Term(block, self)
                self.terms[term.id] = term
            elif block.startswith("[Typedef]"):
                typedef = Typedef(block, self)
                self.typedefs[typedef.id] = typedef
            elif block.startswith("[Instance]"):
                instance = Instance(block, self)
                self.instances[instance.id] = instance
            if progressCallback and i in milestones:
                progressCallback(100.0*i/len(data))
        
        self.aliasMapper = {}
        for id, term in self.terms.items():
            for typeId, parent in term.related:
                self.terms[parent].relatedTo.add((typeId, id))

    def GetDefinedSlimsSubsets(self):
        """ Return a list of defined subsets
        """
        return [line.split()[1] for line in self.header.split("\n") if line.startswith("subsetdef:")]

    def SetSlimsSubset(self, subset):
        """ Set the slims term subset to subset. If subset is a string it must equal one of the defined subsetdef.
        """
        if type(subset) == str:
            self.slimsSubset = [id for id, term in self.terms.items() if subset in getattr(term, "subset", set())]
        else:
            self.slimsSubset = set(subset)
        print self.slimsSubset

    def GetSlimTerms(self, termId):
        """ Return a list of slim terms for termId.
        """
        queue = set([termId])
        visited = set()
        slims = set()
        while queue:
            term = queue.pop()
            visited.add(term)
            if term in self.slimsSubset:
                slims.add(term)
            else:
                queue.update(set(id for typeId, id in self.terms[term].related) - visited)
        return slims

    def ExtractSuperGraph(self, terms):
        """ Return all super terms of terms up to the most general one.
        """
        terms = [terms] if type(terms) == str else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(set(id for typeId, id in self.terms[term].related) - visited)
        return visited

    def ExtractSubGraph(self, terms):
        """ Return all sub terms of terms.
        """
        terms = [terms] if type(terms) == str else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(set(id for typeId, id in self.terms[term].relatedTo) - visited)
        return visited

    def GetTermDepth(self, term, cache_={}):
        """ Return the minimum depth of a term (length of the shortest path to this term from the top level term).
        """
        if term not in cache:
            cache[term] = min([self.GetTermDepth(parent) + 1 for typeId, parent in self.terms[term].related] or [1])
        return cache[term]

    def __getitem__(self, id):
        """ Return object with id (same as ontology.terms[id]
        """
        return self.terms[id]

    def __iter__(self):
        """ Iterate over all ids in ontology
        """
        return iter(self.terms)

    def __len__(self):
        """ Return number of objects in ontology
        """
        return len(self.terms)

    def __contains__(self, id):
        return id in self.terms

    @staticmethod
    def DownloadOntology(file, progressCallback=None):
        tFile = tarfile.open(file, "w:gz") if type(file) == str else file
        tmpDir = os.path.join(orngEnviron.bufferDir, "tmp_go/")
        try:
            os.mkdir(tmpDir)
        except Exception:
            pass
        urlretrieve("http://www.geneontology.org/ontology/gene_ontology_edit.obo", os.path.join(tmpDir, "gene_ontology_edit.obo"), progressCallback and __progressCallbackWrapper(progressCallback))
        tFile.add(os.path.join(tmpDir, "gene_ontology_edit.obo"), "gene_ontology_edit.obo")
        tFile.close()
        os.remove(os.path.join(tmpDir, "gene_ontology_edit.obo"))

_re_obj_name_ = re.compile("([a-zA-z0-9-_]+)")

class AnnotationRecord(object):
    """Holds the data for an annotation record read from the annotation file. Fields can be
    accessed with the names: DB, DB_Object_ID, DB_Object_Symbol, Qualifier, GO_ID, DB_Reference,
    Evidence_code, With_or_From, Aspect, DB_Object_Name, DB_Object_Synonym, DB_Object_Type, taxon,
    Date, Assigned_by (e.g. rec.GO_ID)
    or by supplying the original name of the field (see http://geneontology.org/GO.annotation.shtml#file)
    to the get method (e.g. rec.get("GO ID"))
    The object also provides the folowing data members for quicker access: geneName, GOId, evidence,
    aspect and alias(a list of aliases)
    """
    __slots__ = ["original", "geneName", "GOId", "evidence", "aspect", "alias", "aditionalAliases"]
    def __init__(self, fullText):
        self.original = tuple([t.strip() for t in fullText.split("\t")])
        self.geneName = self.original[2]
        self.GOId = self.original[4]
        self.evidence = self.original[6]
        self.aspect = self.original[8]
        self.alias = self.original[10].split("|")
##        for key, val in zip(annotationFields, self.original):
##            self.__dict__[key] = val

        self.aditionalAliases = []
        if ":" in self.DB_Object_Name:
            self.aditionalAliases = _re_obj_name_.findall(self.DB_Object_Name.split(":")[0])

    def __getattr__(self, name):
        if name in annotationFieldsDict:
            return self.original[annotationFieldsDict[name]]
        else:
            raise AttributeError(name)


class Annotations(object):
    """Annotations object holds the annotations.
    """
    def __init__(self, file=None, ontology=None, progressCallback=None):
        """Initialize the annotations from file by calling ParseFile on it. The ontology must be an instance of Ontology class. The optional progressCallback will be called with a single argument to report on the progress.
        """
        self.file = file
        self.ontology = ontology
        self.allAnnotations = defaultdict(list)
        self.geneAnnotations = defaultdict(list)
        self.termAnnotations = defaultdict(list)
        self.geneNames = set()
        self.geneNamesDict = {}
        self.aliasMapper = {}
        self.additionalAliases = {}
        self.annotations = []
        self.header = ""
        self.geneMapper = None
        if file and os.path.exists(file):
            self.ParseFile(file, progressCallback)
        elif file:
            a = self.Load(file, ontology, progressCallback)
            self.__dict__ = a.__dict__

    def SetOntology(self, ontology):
        self.allAnnotations = defaultdict(list)
        self._ontology = ontology

    def GetOntology(self):
        return self._ontology

    ontology = property(GetOntology, SetOntology, doc="Ontology object for annotations")

    @classmethod
    def Load(cls, org, ontology=None, progressCallback=None):
        """A class method that tries to load the association file for the given organism from default_database_path.
        """
##        files = [name for name in os.listdir(default_database_path) if name.startswith("gene_association") and org in name and not name.endswith(".info")]
        import orngServerFiles
        ## Check if given org is a matching code for GO organism
        files = [file for file in orngServerFiles.listfiles("GO") if file.startswith("gene_association") and org in file.split(".")[:2]]
        name = org
        if not files:
            print >> sys.stderr, "Annotation file not found on local drive\nSearching for annotations on server"
            serverFiles = orngServerFiles.ServerFiles()
            files = [file for file in serverFiles.listfiles("GO") if file.startswith("gene_association") and org in file.split(".")[:2]]
        if not files:
            print >> sys.stderr, "Unable to find annotations for", org, "Matching name against NCBI Taxonomy"
            import obiTaxonomy as tax
            ids = tax.search(org)
            ids = set(ids).intersection(Taxonomy().tax.keys())
            codes = reduce(set.union, [from_taxid(id) for id in ids], set())
            if len(codes) > 1:
                raise tax.MultipleSpeciesException, ", ".join(["%s: %s" % (str(from_taxid(id)), tax.name(id)) for id in ids])
            elif len(codes) == 0:
                raise tax.UnknownSpeciesIdentifier, org
            name, code = tax.name(ids.pop()), codes.pop()
            print >> sys.stderr, "Found annotations for", name, "(%s)" % code
            files = ["gene_association.%s.tar.gz" % code]

        path = os.path.join(orngServerFiles.localpath("GO"), files[0])
        if not os.path.exists(path):
            print >> sys.stderr, "Downloading annotations for", name
            orngServerFiles.download("GO", files[0])
            
        return cls(path, ontology=ontology, progressCallback=progressCallback)
    
    def ParseFile(self, file, progressCallback=None):
        """ Parse and load the annotations from file. Report progress with progressCallback.
        File can be:
            - a tarball containing the association file named gene_association
            - a directory name containing the association file named gene_association
            - a path to the actual association file
            - an open file-like object of the association file
        """
        if type(file) == str:
            if os.path.isfile(file) and tarfile.is_tarfile(file):
                f = tarfile.open(file).extractfile("gene_association")
            elif os.path.isfile(file):
                f = open(file)
            else:
                f = open(os.path.join(file, "gene_association"))
        else:
            f = file
        lines = [line for line in f.read().split("\n") if line.strip()]
        milestones = set(i for i in range(0, len(lines), max(len(lines)/100, 1)))
        for i,line in enumerate(lines):
            if line.startswith("!"):
                self.header = self.header + line + "\n"
                continue
            
            a=AnnotationRecord(line)
            if not a.geneName or not a.GOId or a.Qualifier == "NOT":
                continue
            if a.geneName not in self.geneNames:
                self.geneNames.add(a.geneName)
                self.geneAnnotations[a.geneName].append(a)
                for alias in a.alias:
                    self.aliasMapper[alias] = a.geneName
                for alias in a.aditionalAliases:
                    self.additionalAliases[alias] = a.geneName
                self.aliasMapper[a.geneName] = a.geneName
                self.aliasMapper[a.DB_Object_ID] = a.geneName
                names = [a.DB_Object_ID, a.DB_Object_Symbol]
                names.extend(a.alias)
                for n in names:
                    self.geneNamesDict[n] = names
            else:
                self.geneAnnotations[a.geneName].append(a)
            self.annotations.append(a)
            self.termAnnotations[a.GOId].append(a)
            if progressCallback and i in milestones:
                progressCallback(100.0*i/len(lines))

    def GetGeneNamesTranslator(self, genes):
        def alias(gene):
            return self.aliasMapper.get(gene, self.additionalAliases.get(gene, None))
        return dict([(alias(gene), gene) for gene in genes if alias(gene)])

    def _CollectAnnotations(self, id):
        """ Recursive function collects and caches all annotations for id
        """
        if id not in self.allAnnotations:
            annotations = [self.termAnnotations[id]] ## annotations for this term alone
            for typeId, child in self.ontology[id].relatedTo:
                aa = self._CollectAnnotations(child)
                if type(aa) == set: ## if it was allready reduced in GetAllAnnotations
                    annotations.append(aa)
                else:
                    annotations.extend(aa)
            self.allAnnotations[id] = annotations
        return self.allAnnotations[id]

    def GetAllAnnotations(self, id):
        """ Return a set of all annotations for this and all subterms.
        """
        if id not in self.allAnnotations or type(self.allAnnotations[id]) == list:
            annot_set = set()
            for annots in self._CollectAnnotations(id):
                annot_set.update(annots)
            self.allAnnotations[id] = annot_set
        return self.allAnnotations[id]


    def GetAllGenes(self, id, evidenceCodes = None):
        """ Return a list of genes annotated by specified evidence codes to this and all subterms."
        """
        evidenceCodes = set(evidenceCodes or evidenceDict.keys())
        annotations = self.GetAllAnnotations(id)
        return list(set([ann.geneName for ann in annotations if ann.Evidence_code in evidenceCodes]))

    def GetEnrichedTerms(self, genes, reference=None, evidenceCodes=None, slimsOnly=False, aspect="P", prob=obiProb.Binomial(), progressCallback=None):
        """ Return a dictionary of enriched terms, with tuples of (list_of_genes, p_value, reference_count) for items and term ids as keys.
        """
        revGenesDict = self.GetGeneNamesTranslator(genes)
        genes = set(revGenesDict.keys())
        reference = set(reference) if reference else self.geneNames
        evidenceCodes = set(evidenceCodes or evidenceDict.keys())
        annotations = [ann for gene in genes for ann in self.geneAnnotations[gene] if ann.Evidence_code in evidenceCodes and ann.Aspect == aspect]
        refAnnotations = set([ann for gene in reference for ann in self.geneAnnotations[gene] if ann.Evidence_code in evidenceCodes and ann.Aspect == aspect])
        annotationsDict = defaultdict(set)
        for ann in annotations:
            annotationsDict[ann.GO_ID].add(ann)
            
        terms = self.ontology.ExtractSuperGraph(annotationsDict.keys())
        res = {}
        milestones = set(range(0, len(terms), max(len(terms)/100, 1)))
        for i, term in enumerate(terms):
            if slimsOnly and term not in self.ontology.slimsSubset:
                continue
            allAnnotations = self.GetAllAnnotations(term)
            allAnnotations.intersection_update(refAnnotations)
            allAnnotatedGenes = set([ann.geneName for ann in allAnnotations])
            if len(genes) > len(allAnnotatedGenes): 
                mappedGenes = genes.intersection(allAnnotatedGenes)
            else:
                mappedGenes = allAnnotatedGenes.intersection(genes)
            if len(reference) > len(allAnnotatedGenes):
                mappedReferenceGenes = reference.intersection(allAnnotatedGenes)
            else:
                mappedReferenceGenes = allAnnotatedGenes.intersection(reference)
            res[term] = ([revGenesDict[g] for g in mappedGenes], prob.p_value(len(mappedGenes), len(reference), len(mappedReferenceGenes), len(genes)), len(mappedReferenceGenes))
            if progressCallback and i in milestones:
                progressCallback(100.0 * i / len(terms))
        return res

    def GetAnnotatedTerms(self, genes, directAnnotationOnly=False, evidenceCodes=None, progressCallback=None):
        """ Return all terms that are annotated by genes with evidenceCodes.
        """
        genes = [genes] if type(genes) == str else genes
        revGenesDict = self.GetGeneNamesTranslator(genes)
        genes = set(revGenesDict.keys())
        evidenceCodes = set(evidenceCodes or evidenceDict.keys())
        annotations = [ann for gene in genes for ann in self.geneAnnotations[gene] if ann.Evidence_code in evidenceCodes]
        dd = defaultdict(set)
        for ann in annotations:
            dd[ann.GO_ID].add(ann.geneName)
        if not directAnnotationOnly:
            terms = self.ontology.ExtractSuperGraph(dd.keys())
            for i, term in enumerate(terms):
                termAnnots = self.GetAllAnnotations(term)
                termAnnots.intersection_update(annotations)
                dd[term].update([revGenesDict.get(ann.geneName, ann.geneName) for ann in termAnnots])
        return dict(dd)

    def __iter__(self):
        """ Iterate over all AnnotationRecord objects in annotations
        """
        return iter(self.annotations)

    def __len__(self):
        """ Return the number of annotations
        """
        return len(self.annotations)

    def __getitem__(self, index):
        """ Return the i-th annotation record
        """
        return self.annotations[index]
    
    @staticmethod
    def DownloadAnnotations(org, file, progressCallback=None):
        tFile = tarfile.open(file, "w:gz") if type(file) == str else file
        tmpDir = os.path.join(orngEnviron.bufferDir, "tmp_go/")
        try:
            os.mkdir(tmpDir)
        except Exception:
            pass
        fileName = "gene_association." + org + ".gz"
        urlretrieve("http://www.geneontology.org/gene-associations/" + fileName, os.path.join(tmpDir, fileName), progressCallback and __progressCallbackWraper(progressCallback))
        gzFile = GzipFile(os.path.join(tmpDir, fileName), "r")
        file = open(os.path.join(tmpDir, "gene_association." + org), "w")
        file.writelines(gzFile.readlines())
        file.flush()
        file.close()
##        tFile = tarfile.open(os.path.join(tmpDir, "gene_association." + org + ".tar.gz"), "w:gz")
        tFile.add(os.path.join(tmpDir, "gene_association." + org), "gene_association")
        annotation = Annotations(os.path.join(tmpDir, "gene_association." + org), progressCallback=progressCallback)
        cPickle.dump(annotation.geneNames, open(os.path.join(tmpDir, "gene_names.pickle"), "wb"))
        tFile.add(os.path.join(tmpDir, "gene_names.pickle"), "gene_names.pickle")
        tFile.close()
        os.remove(os.path.join(tmpDir, "gene_association." + org))
        os.remove(os.path.join(tmpDir, "gene_names.pickle"))

class Taxonomy(object):
    """Maps NCBI taxonomy ids to coresponding GO organism codes
    """
    __shared_state = {"tax": None}
    def __init__(self):
        self.__dict__ = self.__shared_state
        if not self.tax:
            import orngServerFiles
            path = orngServerFiles.localpath("GO", "taxonomy.pickle")
            if os.path.isfile(path):
                self.tax = cPickle.load(open(path, "rb"))
            else:
                orngServerFiles.download("GO", "taxonomy.pickle")
                self.tax = cPickle.load(open(path, "rb"))
                
    def __getitem__(self, key):
        return list(self.tax[key])
    
def from_taxid(id):
    """ Return a set of GO organism codes that correspond to NCBI taxonomy id
    """
    return set(Taxonomy()[id])

def to_taxid(db_code):
    """ Return a set of NCBI taxonomy ids from db_code GO organism annotations
    """
    r = [key for key, val in Taxonomy().tax.items() if db_code in val]
    return set(r)
    

class __progressCallbackWrapper:
    def __init__(self, callback):
        self.callback = callback
    def __call__(self, bCount, bSize, fSize):
        fSize = 10000000 if fSize == -1 else fSize
        self.callback(100*bCount*bSize/fSize)
        
from obiGenomicsUpdate import Update as UpdateBase

import urllib2

class Update(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path or getDataDir(), progressCallback)
    def CheckModified(self, addr, date=None):
        return date > self.GetLastModified(addr) if date else True
        
    def CheckModifiedOrg(self, org):
        return self.CheckModified("http://www.geneontology.org/gene-associations/gene_association." + org + ".gz", self.LastModifiedOrg(org))
    
    def LastModifiedOrg(self, org):
        return self.shelve.get((Update.UpdateAnnotation, (org,)), None)

    def GetLastModified(self, addr):
        stream = urllib2.urlopen(addr)
        return datetime.strptime(stream.headers.get("Last-Modified"), "%a, %d %b %Y %H:%M:%S %Z")
##        return stream.headers.get("Last-Modified")

    def GetAvailableOrganisms(self):
        source = urllib2.urlopen("http://www.geneontology.org/gene-associations/").read()
        return [s.split(".")[1] for s in sorted(set(re.findall("gene_association\.[a-zA-z0-9_]+?\.gz", source)))]

    def GetDownloadedOrganisms(self):
        return [name.split(".")[1] for name in os.listdir(self.local_database_path) if name.startswith("gene_association")]

    def IsUpdatable(self, func, args):
        if func == Update.UpdateOntology:
            return self.CheckModified("http://www.geneontology.org/ontology/gene_ontology.obo", self.shelve.get((Update.UpdateOntology, ()), None))
        elif func == Update.UpdateAnnotation:
            return self.CheckModifiedOrg(args[0])
            
    def GetDownloadable(self):
        orgs = set(self.GetAvailableOrganisms()) - set(self.GetDownloadedOrganisms())
        ret = []
        if (Update.UpdateOntology, ()) not in self.shelve:
            ret.append((Update.UpdateOntology, ()))
        if orgs:
            ret.extend([(Update.UpdateAnnotation, (org,)) for org in orgs])
        return ret

    def UpdateOntology(self):
        Ontology.DownloadOntology(os.path.join(self.local_database_path, "gene_ontology_edit.obo.tar.gz"), self.progressCallback)
        self._update(Update.UpdateOntology, (), self.GetLastModified("http://www.geneontology.org/ontology/gene_ontology.obo"))

    def UpdateAnnotation(self, org):
        Annotations.DownloadAnnotations(org, os.path.join(self.local_database_path, "gene_association." + org + ".tar.gz"), self.progressCallback)
        self._update(Update.UpdateAnnotation, (org,), self.GetLastModified("http://www.geneontology.org/gene-associations/gene_association." + org + ".gz"))
        
    def UpdateTaxonomy(self, org):
        exclude = ["goa_uniprot", "goa_pdb", "GeneDB_tsetse", "reactome", "goa_zebrafish", "goa_rat", "goa_mouse"]

        orgs = self.GetAvailableOrganisms()
        tax = defaultdict(set)

        for org in orgs:
            if org in exclude:
                continue
            try:
                a = obiGO.Annotations(os.path.join(self.local_database_path, "gene_association." + org + ".tar.gz"))
                taxons = set(ann.taxon for ann in a.annotations)
                for taxId in [t.split(":")[-1] for t in taxons if "|" not in t]: ## exclude taxons with cardinality 2
                    tax[taxId].add(org)
            except Exception, ex:
                print ex
                
        cPickle.dump(dict(tax), open(os.path.join(path, "taxonomy.pickle"), "wb"))
            

def _test1():
##    Ontology.DownloadOntology("ontology_arch.tar.gz")
##    Annotations.DownloadAnnotations("sgd", "annotations_arch.tar.gz")
    def _print(f):
        print f
    o = Ontology("ontology_arch.tar.gz")
    a = Annotations("annotations_arch.tar.gz", ontology=o)
    import go
    loadAnnotation("sgd")
    loadGO()
##    print a.GetAllGenes("GO:0008150")
    import profile
    a.GetEnrichedTerms(sorted(a.geneNames)[:100])#, progressCallback=_print)
##    profile.runctx("a.GetEnrichedTerms(sorted(a.geneNames)[:100])", {"a":a}, {})
    a.GetEnrichedTerms(sorted(a.geneNames)[:100])#, progressCallback=_print)
    d1 = a.GetEnrichedTerms(sorted(a.geneNames)[:1000])#, progressCallback=_print)
    d2 = GOTermFinder(sorted(a.geneNames)[:1000])
    print set(d2.keys()) - set(d1.keys())
    print set(d1.keys()) - set(d2.keys())
##    print a.GetEnrichedTerms(sorted(a.geneNames)[:100])#, progressCallback=_print)
    
if __name__ == "__main__":
    _test1()
