"""
`obiGO` is a Gene Ontology (GO) Handling Library.

"""

from __future__ import absolute_import
import os
import sys
import tarfile
import gzip
import re
import cPickle
import shutil
import urllib2
import warnings
import copy

from urllib import urlretrieve
from gzip import GzipFile
from datetime import datetime
from collections import defaultdict

from Orange.orng import orngEnviron
from Orange.orng import orngServerFiles
from Orange.orng import orngMisc

from . import obiProb, obiGene, obiTaxonomy

default_database_path = os.path.join(orngServerFiles.localpath(), "GO")


_CVS_REVISION_RE = re.compile(r"^(rev)?\d+\.\d+$")

evidenceTypes = {
# Experimental
    'EXP': 'Inferred from Experiment',
    'IDA': 'Inferred from Direct Assay',
    'IPI': 'Inferred from Physical Interaction',  # [with <database:protein_name>]',
    'IMP': 'Inferred from Mutant Phenotype',
    'IGI': 'Inferred from Genetic Interaction',  # [with <database:gene_symbol[allele_symbol]>]',
    'IEP': 'Inferred from Expression Pattern',
# Computational Analysis Evidence Codes
    'ISS': 'Inferred from Sequence Similarity',  # [with <database:sequence_id>] ',
    'ISA': 'Inferred from Sequence Alignment',
    'ISO': 'Inferred from Sequence Orthology',
    'ISM': 'Inferred from Sequence Model',
    'IGC': 'Inferred from Genomic Context',
    'RCA': 'Inferred from Reviewed Computational Analysis',
# Author Statement Evidence Codes
    'TAS': 'Traceable author statement',
    'NAS': 'Non-traceable author statement',
# Curatorial Statement Evidence Codes
    'IC': 'Inferred by curator',
    'ND': 'No biological data available',
# Computationally-assigned Evidence Codes
    'IEA': 'Inferred from electronic annotation',  # [to <database:id>]',
# Obsolete Evidence Codes
    'NR': 'Not Recorded(Obsolete)'
}


evidenceDict = defaultdict(int, [(e, 2 ** i) for i, e in
                                 enumerate(evidenceTypes.keys())])

evidenceTypesOrdered = [
'EXP',
'IDA',
'IPI',
'IMP',
'IGI',
'IEP',
# Computational Analysis Evidence Codes
'ISS',
'ISA',
'ISO',
'ISM',
'IGC',
'RCA',
# Author Statement Evidence Codes
'TAS',
'NAS',
# Curatorial Statement Evidence Codes
'IC',
'ND',
# Computationally-assigned Evidence Codes
'IEA',
# Obsolete Evidence Codes
'NR'
]

multiplicitySet = set(
    ["alt_id", "is_a", "subset", "synonym", "related_synonym",
     "exact_synonym", "broad_synonym", "narrow_synonym",
     "xref_analog", "xref_unknown", "relationship"])

multipleTagSet = multiplicitySet

annotationFields = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
    "DB_Reference", "Evidence_Code", "With_From", "Aspect", "DB_Object_Name",
    "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_by"]

annotationFieldsDict = {
    "DB": 0,
    "DB_Object_ID": 1,
    "DB_Object_Symbol": 2,
    "Qualifier": 3,
    "GO_ID": 4,
    "GO ID": 4,
    "GOID": 4,
    "DB_Reference": 5,
    "DB:Reference": 5,
    "Evidence_Code": 6,
    "Evidence Code": 6,
    "Evidence_code": 6,  # compatibility with  older revisions
    "With_or_From": 7,
    "With (or) From": 7,
    "Aspect": 8,
    "DB_Object_Name": 9,
    "DB_Object_Synonym": 10,
    "DB_Object_Type": 11,
    "taxon": 12,
    "Date": 13,
    "Assigned_by": 14
}

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
    """Represents a generic OBO object (e.g. Term, Typedef, Instance, ...)
    Example:

    >>> OBOObject(r"[Term]\nid: FOO:001\nname: bar", ontology)

    """
    _INTERN_TAGS = ["id", "name", "namespace", "alt_id", "is_a"]

    def __init__(self, stanza=None, ontology=None):
        self.ontology = ontology
        self._lines = []
        self.values = {}
        self.related = set()
        self.relatedTo = set()
        if stanza:
            self.ParseStanza(stanza)

    def ParseStanza(self, stanza):
        intern_tags = set(self._INTERN_TAGS)
        for line in stanza.splitlines():
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
            tag = intern(tag)
            value = value.strip()
            comment = comment.strip()
            if tag in intern_tags:
                value, comment = intern(value), intern(comment)
            self._lines.append((tag, value, modifiers, comment))
            if tag in multipleTagSet:
                self.values.setdefault(tag, []).append(value)
            else:
                self.values[tag] = value
        self.related = set(self.GetRelatedObjects())
        self.__dict__.update(self.values)
        if "def" in self.__dict__:
            self.__dict__["def_"] = self.__dict__["def"]

    def GetRelatedObjects(self):
        """Return a list of tuple pairs where the first element is relationship
        typeId and the second id of object to whom the relationship applies to.

        """
        # TODO: add other defined Typedef ids
        typeIds = [intern("is_a")]
        result = [(typeId, id) for typeId in typeIds
                  for id in self.values.get(typeId, [])]
        result = result + [tuple(map(intern, r.split(None, 1)))
                           for r in self.values.get("relationship", [])]
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
    """ Ontology is the main class representing a gene ontology.

    Example::
        >>> ontology = Ontology("my_gene_ontology.obo")

    """
    version = 1

    def __init__(self, file=None, progressCallback=None, rev=None):
        """ Initialize the ontology from file (if `None` the default gene
        ontology will be loaded). The optional `progressCallback` will be
        called with a single argument to report on the progress.

        """
        self.terms = {}
        self.typedefs = {}
        self.instances = {}
        self.slimsSubset = set()

        if file is not None:
            self.ParseFile(file, progressCallback)
        elif rev is not None:
            if not _CVS_REVISION_RE.match(rev):
                raise ValueError("Invalid revision format.")
            if rev.startswith("rev"):
                rev = rev[3:]
            pc = lambda v: progressCallback(v / 2.0) \
                 if progressCallback else None
            filename = os.path.join(default_database_path,
                                    "gene_ontology_edit@rev%s.obo" % rev)
            if not os.path.exists(filename):
                self.DownloadOntologyAtRev(rev, filename, pc)
            self.ParseFile(filename, lambda v: progressCallback(v / 2.0 + 50) \
                                     if progressCallback else None)
        else:
            # Load the default ontology file.
            fool = self.Load(progressCallback)
            # A fool and his attributes are soon parted
            self.__dict__ = fool.__dict__

    @classmethod
    def Load(cls, progressCallback=None):
        """ A class method that tries to load the ontology file from
        default_database_path. It looks for a filename starting with
        'gene_ontology'. If not found it will download it.

        """
        filename = os.path.join(default_database_path,
                                "gene_ontology_edit.obo.tar.gz")
        if not os.path.isfile(filename) and not os.path.isdir(filename):
            orngServerFiles.download("GO", "gene_ontology_edit.obo.tar.gz")

        return cls(filename, progressCallback=progressCallback)

    def ParseFile(self, file, progressCallback=None):
        """ Parse the file. file can be a filename string or an open filelike
        object. The optional progressCallback will be called with a single
        argument to report on the progress.
        """
        if isinstance(file, basestring):
            if os.path.isfile(file) and tarfile.is_tarfile(file):
                f = tarfile.open(file).extractfile("gene_ontology_edit.obo")
            elif os.path.isfile(file):
                f = open(file)
            elif os.path.isdir(file):
                f = open(os.path.join(file, "gene_ontology_edit.obo"))
            else:
                raise ValueError("Cannot open %r for parsing" % file)
        else:
            f = file

        data = f.readlines()
        data = "".join([line for line in data if not line.startswith("!")])
        self.header = data[: data.index("[Term]")]
        c = re.compile("\[.+?\].*?\n\n", re.DOTALL)
        data = c.findall(data)

        milestones = orngMisc.progressBarMilestones(len(data), 90)
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
                progressCallback(90.0 * i / len(data))

        self.aliasMapper = {}
        self.reverseAliasMapper = defaultdict(set)
        milestones = orngMisc.progressBarMilestones(len(self.terms), 10)
        for i, (id, term) in enumerate(self.terms.iteritems()):
            for typeId, parent in term.related:
                self.terms[parent].relatedTo.add((typeId, id))
            try:
                self.aliasMapper.update([(alt_id, id)
                                         for alt_id in term.alt_id])
                self.reverseAliasMapper[id].union_update(term.alt_id)
            except AttributeError:
                pass
            if progressCallback and i in milestones:
                progressCallback(90.0 + 10.0 * i / len(self.terms))

    def GetDefinedSlimsSubsets(self):
        """ Return a list of defined subsets.
        """
        return [line.split()[1] for line in self.header.splitlines()
                if line.startswith("subsetdef:")]

    def SetSlimsSubset(self, subset):
        """ Set the slims term subset to subset. If subset is a string it
        must equal one of the defined subsetdef.
        """
        if type(subset) == str:
            self.slimsSubset = [id for id, term in self.terms.items()
                                if subset in getattr(term, "subset", set())]
        else:
            self.slimsSubset = set(subset)

    def GetSlimsSubset(self, subset):
        return [id for id, term in self.terms.items()
                if subset in getattr(term, "subset", set())]

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
                queue.update(set(id for typeId, id in self[term].related) -
                             visited)
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
            queue.update(set(id for typeId, id in self[term].related) -
                         visited)
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
            queue.update(set(id for typeId, id in self[term].relatedTo) -
                         visited)
        return visited

    def GetTermDepth(self, term, cache_={}):
        """ Return the minimum depth of a term (length of the shortest
        path to this term from the top level term).
        """
        if term not in cache_:
            cache_[term] = min([self.GetTermDepth(parent) + 1
                                for typeId, parent in self[term].related] or
                               [1])
        return cache_[term]

    def __getitem__(self, id):
        """
        Return a object with id.
        """
        if id in self.terms:
            return self.terms[id]
        elif id in self.aliasMapper:
            return self.terms[self.aliasMapper[id]]
        else:
            raise KeyError(id)

    def __iter__(self):
        """
        Iterate over all ids in ontology.
        """
        return iter(self.terms)

    def __len__(self):
        """
        Return number of objects in ontology.
        """
        return len(self.terms)

    def __contains__(self, id):
        """
        Return `True` if a term with `id` is present in the ontology.
        """
        return id in self.terms or id in self.aliasMapper

    @staticmethod
    def DownloadOntology(file, progressCallback=None):
        tFile = tarfile.open(file, "w:gz") if isinstance(file, basestring) \
                else file
        tmpDir = os.path.join(orngEnviron.bufferDir, "tmp_go/")
        try:
            os.mkdir(tmpDir)
        except Exception:
            pass
        urlretrieve("http://www.geneontology.org/ontology/gene_ontology_edit.obo",
                    os.path.join(tmpDir, "gene_ontology_edit.obo"),
                    progressCallback and __progressCallbackWrapper(progressCallback))
        tFile.add(os.path.join(tmpDir, "gene_ontology_edit.obo"),
                  "gene_ontology_edit.obo")
        tFile.close()
        os.remove(os.path.join(tmpDir, "gene_ontology_edit.obo"))

    @staticmethod
    def DownloadOntologyAtRev(rev, filename=None, progressCallback=None):
        url = "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/ontology/gene_ontology_edit.obo?rev=%s" % rev
        url += ";content-type=text%2Fplain"
        if filename is None:
            filename = os.path.join(default_database_path,
                                    "gene_ontology_edit@rev%s.obo" % rev)
        r = urllib2.urlopen(url)

        with open(filename + ".part", "wb") as f:
            shutil.copyfileobj(r, f)

        os.rename(filename + ".part", filename)


_re_obj_name_ = re.compile("([a-zA-z0-9-_]+)")


class AnnotationRecord(object):
    """Holds the data for an annotation record read from the annotation file.
    Fields can be accessed with the names: DB, DB_Object_ID, DB_Object_Symbol,
    Qualifier, GO_ID, DB_Reference, Evidence_code, With_or_From, Aspect,
    DB_Object_Name, DB_Object_Synonym, DB_Object_Type, taxon, Date,
    Assigned_by (e.g. rec.GO_ID) or by supplying the original name of the
    field (see http://geneontology.org/GO.annotation.shtml#file) to the get
    method (e.g. rec.get("GO ID")) The object also provides the following
    data members for quicker access: geneName, GOId, evidence, aspect and
    alias(a list of aliases)

    """
    __slots__ = annotationFields + ["geneName", "GOId", "evidence",
                                    "aspect", "alias", "additionalAliases"]

    def __init__(self, fullText):
        """\
        :param fulText: A single line from the annotation file.

        """
        for slot, val in zip(annotationFields, fullText.split("\t")):
            setattr(self, slot, intern(val))
        self.geneName = self.DB_Object_Symbol
        self.GOId = self.GO_ID
        self.evidence = self.Evidence_Code
        self.aspect = self.Aspect
        self.alias = list(map(intern, self.DB_Object_Synonym.split("|")))

        self.additionalAliases = []
        if ":" in self.DB_Object_Name:
            self.additionalAliases = []  # _re_obj_name_.findall(self.DB_Object_Name.split(":")[0])

    def __getattr__(self, name):
        if name in annotationFieldsDict:
            return getattr(self, self.__slots__[annotationFieldsDict[name]])
        else:
            raise AttributeError(name)


class Annotations(object):
    """Annotations object holds the annotations.
    """
    version = 2

    def __init__(self, file=None, ontology=None, genematcher=None,
                 progressCallback=None, rev=None):
        """Initialize the annotations from file by calling ParseFile on it.
        The ontology must be an instance of Ontology class.
        The optional progressCallback will be called with a single argument
        to report on the progress.

        """
        self.ontology = ontology
        self.allAnnotations = defaultdict(list)
        self.geneAnnotations = defaultdict(list)
        self.termAnnotations = defaultdict(list)
        self._geneNames = None
        self._geneNamesDict = None
        self._aliasMapper = None
        self.additionalAliases = {}
        self.annotations = []
        self.header = ""
        self.genematcher = genematcher
        self.taxid = None

        if type(file) in [list, set, dict, Annotations]:
            for ann in file:
                self.AddAnnotation(ann)
            if type(file, Annotations):
                self.taxid = file.taxid

        elif isinstance(file, basestring) and os.path.exists(file):
            self.ParseFile(file, progressCallback)
            try:
                self.taxid = to_taxid(os.path.basename(file).split(".")[1]).pop()
            except IOError:
                pass

        elif isinstance(file, basestring):
            # Assuming organism code/name
            if rev is not None:
                if not _CVS_REVISION_RE.match(rev):
                    raise ValueError("Invalid revision format")

                if rev.startswith("rev"):
                    rev = rev[3:]
                code = self.organism_name_search(file)
                filename = os.path.join(default_database_path,
                                        "gene_association.%s@rev%s.tar.gz" %
                                        (code, rev))

                if not os.path.exists(filename):
                    self.DownloadAnnotationsAtRev(
                        code, rev, filename, progressCallback)

                self.ParseFile(filename, progressCallback)
                self.taxid = to_taxid(code).pop()
            else:
                a = self.Load(file, ontology, genematcher, progressCallback)
                self.__dict__ = a.__dict__
                self.taxid = to_taxid(organism_name_search(file)).pop()
        elif file is not None:
            self.ParseFile(file, progressCallback)

        if not self.genematcher and self.taxid:
            matchers = [obiGene.GMGO(self.taxid)]
            if self.taxid == "352472":
                matchers.extend(
                    [obiGene.GMDicty(),
                     [obiGene.GMGO(self.taxid), obiGene.GMDicty()]]
                )

            self.genematcher = obiGene.matcher(matchers)

        if self.genematcher:
            self.genematcher.set_targets(self.geneNames)

    @classmethod
    def organism_name_search(cls, org):
        ids = to_taxid(org)
        from . import obiTaxonomy as tax
        if not ids:
            ids = [org] if org in Taxonomy().common_org_map.keys() + \
                  Taxonomy().code_map.keys() else []
        if not ids:
            ids = tax.to_taxid(org, mapTo=Taxonomy().keys())
        if not ids:
            ids = tax.search(org, exact=True)
            ids = set(ids).intersection(Taxonomy().keys())
        if not ids:
            ids = tax.search(org)
            ids = set(ids).intersection(Taxonomy().keys())
        codes = set([from_taxid(id) for id in ids])
        if len(codes) > 1:
            raise tax.MultipleSpeciesException(
                ", ".join(["%s: %s" % (str(from_taxid(id)), tax.name(id))
                           for id in ids]))
        elif len(codes) == 0:
            raise tax.UnknownSpeciesIdentifier(org)
        return codes.pop()

    @classmethod
    def organism_version(cls, name):
        name = organism_name_search(name)
        orngServerFiles.localpath_download(
            "GO", "gene_association.%s.tar.gz" % name)
        return ("v%i." % cls.version) + orngServerFiles.info("GO",
                        "gene_association.%s.tar.gz" % name)["datetime"]

    def SetOntology(self, ontology):
        """ Set the ontology to use in the annotations mapping.
        """
        self.allAnnotations = defaultdict(list)
        self._ontology = ontology

    def GetOntology(self):
        return self._ontology

    ontology = property(GetOntology, SetOntology,
                        doc="Ontology object for annotations")

    @classmethod
    def Load(cls, org, ontology=None, genematcher=None, progressCallback=None):
        """A class method that tries to load the association file for the
        given organism from default_database_path.
        """
        code = organism_name_search(org)

        print "CODE: %s" % code
        print "ORG: %s" % org

        file = "gene_association.%s.tar.gz" % code

        path = os.path.join(orngServerFiles.localpath("GO"), file)

        if not os.path.exists(path):
            sf = orngServerFiles.ServerFiles()
            available = sf.listfiles("GO")
            if file not in available:
                from . import obiKEGG
                raise obiTaxonomy.UnknownSpeciesIdentifier(org + str(code))
            orngServerFiles.download("GO", file)

        return cls(path, ontology=ontology, genematcher=genematcher,
                   progressCallback=progressCallback)

    def ParseFile(self, file, progressCallback=None):
        """ Parse and load the annotations from file. Report progress
        with progressCallback.
        File can be:
            - a tarball containing the association file named gene_association
            - a directory name containing the association file named gene_association
            - a path to the actual association file
            - an open file-like object of the association file
        """
        if isinstance(file, basestring):
            if os.path.isfile(file) and tarfile.is_tarfile(file):
                f = tarfile.open(file).extractfile("gene_association")
            elif os.path.isfile(file) and file.endswith(".gz"):
                f = gzip.open(file)
            elif os.path.isfile(file):
                f = open(file)
            elif os.path.isdir(file):
                f = open(os.path.join(file, "gene_association"))
            else:
                raise ValueError("Cannot open %r for parsing." % file)
        else:
            f = file
        lines = [line for line in f.read().splitlines() if line.strip()]

        milestones = orngMisc.progressBarMilestones(len(lines), 100)
        for i, line in enumerate(lines):
            if line.startswith("!"):
                self.header = self.header + line + "\n"
                continue

            a = AnnotationRecord(line)
            self.AddAnnotation(a)
#            self.annotations.append(a)
            if progressCallback and i in milestones:
                progressCallback(100.0 * i / len(lines))

    def AddAnnotation(self, a):
        """ Add a single `AnotationRecord` instance to this `Annotations`
        object.
        """
        if not isinstance(a, AnnotationRecord):
            a = AnnotationRecord(a)
        if not a.geneName or not a.GOId or a.Qualifier == "NOT":
            return

        self.geneAnnotations[a.geneName].append(a)
        self.annotations.append(a)
        self.termAnnotations[a.GOId].append(a)
        self.allAnnotations = defaultdict(list)

        self._geneNames = None
        self._geneNamesDict = None
        self._aliasMapper = None

    @property
    def geneNamesDict(self):
        if getattr(self, "_geneNamesDict", None) is None:
            self._geneNamesDict = defaultdict(set)
            for alias, name in self.aliasMapper.iteritems():
                self._geneNamesDict[name].add(alias)
        return self._geneNamesDict

    @property
    def geneNames(self):
        if getattr(self, "_geneNames", None) is None:
            self._geneNames = set([ann.geneName for ann in self.annotations])
        return self._geneNames

    @property
    def aliasMapper(self):
        if getattr(self, "_aliasMapper", None) is None:
            self._aliasMapper = {}
            for ann in self.annotations:
                self._aliasMapper.update([(alias, ann.geneName)
                                          for alias in ann.alias +
                                           [ann.geneName, ann.DB_Object_ID]])
        return self._aliasMapper

    def GetGeneNamesTranslator(self, genes):
        """ Return a dictionary mapping canonical names (DB_Object_Symbol)
        to `genes`.

        """
        def alias(gene):
            if self.genematcher:
                return self.genematcher.umatch(gene)
            else:
                return gene if gene in self.geneNames else \
                        self.aliasMapper.get(gene,
                             self.additionalAliases.get(gene, None))

        return dict([(alias(gene), gene) for gene in genes if alias(gene)])

    def _CollectAnnotations(self, id, visited):
        """ Recursive function collects and caches all annotations for id
        """
        if id not in self.allAnnotations and id not in visited:
            if id in self.ontology.reverseAliasMapper:
                annotations = [self.termAnnotations.get(alt_id, [])
                               for alt_id in
                               self.ontology.reverseAliasMapper[id]] + \
                              [self.termAnnotations[id]]
            else:
                ## annotations for this term alone
                annotations = [self.termAnnotations[id]]
            visited.add(id)
            for typeId, child in self.ontology[id].relatedTo:
                aa = self._CollectAnnotations(child, visited)
                if type(aa) == set:
                    ## if it was already reduced in GetAllAnnotations
                    annotations.append(aa)
                else:
                    annotations.extend(aa)
            self.allAnnotations[id] = annotations
        return self.allAnnotations[id]

    def GetAllAnnotations(self, id):
        """ Return a set of all annotations (instances if `AnnotationRectord`)
        for GO term `id` and all it's subterms.

        :param id: GO term id
        :type id: str

        """
        visited = set()
        id = self.ontology.aliasMapper.get(id, id)
        if id not in self.allAnnotations or \
                type(self.allAnnotations[id]) == list:
            annot_set = set()
            for annots in self._CollectAnnotations(id, set()):
                annot_set.update(annots)
            self.allAnnotations[id] = annot_set
        return self.allAnnotations[id]

    def GetAllGenes(self, id, evidenceCodes=None):
        """ Return a list of genes annotated by specified `evidenceCodes`
        to GO term 'id' and all it's subterms."

        :param id: GO term id
        :type id: str

        :param evidneceCodes: List of evidence codes to consider when
                              matching annotations to terms.
        :type evidenceCodes: list-of-strings
        """
        evidenceCodes = set(evidenceCodes or evidenceDict.keys())
        annotations = self.GetAllAnnotations(id)
        return list(set([ann.geneName for ann in annotations
                         if ann.Evidence_Code in evidenceCodes]))

    def GetEnrichedTerms(self, genes, reference=None, evidenceCodes=None,
                         slimsOnly=False, aspect=None, prob=obiProb.Binomial(),
                         useFDR=True, progressCallback=None):
        """ Return a dictionary of enriched terms, with tuples of
        (list_of_genes, p_value, reference_count) for items and term
        ids as keys. P-Values are FDR adjusted if useFDR is True (default).

        :param genes: List of genes
        :param reference: list of genes (if None all genes included in the
                          annotations will be used).
        :param evidenceCodes: List of evidence codes to consider.
        :param slimsOnly: If `True` return only slim terms
        :param aspect: Which aspects to use. Use all by default. "P", "F", "C"
            or a set containing these elements.
        """
        revGenesDict = self.GetGeneNamesTranslator(genes)
        genes = set(revGenesDict.keys())
        if reference:
            refGenesDict = self.GetGeneNamesTranslator(reference)
            reference = set(refGenesDict.keys())
        else:
            reference = self.geneNames

        if aspect == None:
            aspects_set = set(["P", "C", "F"])
        elif isinstance(aspect, basestring):
            aspects_set = set([aspect])
        else:
            aspects_set = aspect

        evidenceCodes = set(evidenceCodes or evidenceDict.keys())
        annotations = [ann
                       for gene in genes for ann in self.geneAnnotations[gene]
                       if ann.Evidence_Code in evidenceCodes and
                       ann.Aspect in aspects_set]

        refAnnotations = set(
            [ann
             for gene in reference for ann in self.geneAnnotations[gene]
             if ann.Evidence_Code in evidenceCodes and
             ann.Aspect in aspects_set]
        )

        annotationsDict = defaultdict(set)
        for ann in annotations:
            annotationsDict[ann.GO_ID].add(ann)

        if slimsOnly and not self.ontology.slimsSubset:
            warnings.warn("Unspecified slims subset in the ontology! "
                          "Using 'goslim_generic' subset", UserWarning)
            self.ontology.SetSlimsSubset("goslim_generic")

        terms = annotationsDict.keys()
        filteredTerms = [term for term in terms if term in self.ontology]

        if len(terms) != len(filteredTerms):
            termDiff = set(terms) - set(filteredTerms)
            warnings.warn("%s terms in the annotations were not found in the "
                          "ontology." % ",".join(map(repr, termDiff)),
                          UserWarning)

        terms = self.ontology.ExtractSuperGraph(filteredTerms)
        res = {}

        milestones = orngMisc.progressBarMilestones(len(terms), 100)
        for i, term in enumerate(terms):
            if slimsOnly and term not in self.ontology.slimsSubset:
                continue
            allAnnotations = self.GetAllAnnotations(term).intersection(refAnnotations)
##            allAnnotations.intersection_update(refAnnotations)
            allAnnotatedGenes = set([ann.geneName for ann in allAnnotations])
            mappedGenes = genes.intersection(allAnnotatedGenes)
##            if not mappedGenes:
##                print >> sys.stderr, term, sorted(genes)
##                print >> sys.stderr, sorted(allAnnotatedGenes)
##                return
            
            if len(reference) > len(allAnnotatedGenes):
                mappedReferenceGenes = reference.intersection(allAnnotatedGenes)
            else:
                mappedReferenceGenes = allAnnotatedGenes.intersection(reference)
            res[term] = ([revGenesDict[g] for g in mappedGenes],
                         prob.p_value(len(mappedGenes), len(reference),
                                      len(mappedReferenceGenes), len(genes)),
                         len(mappedReferenceGenes))
            if progressCallback and i in milestones:
                progressCallback(100.0 * i / len(terms))
        if useFDR:
            res = sorted(res.items(), key=lambda (_1, (_2, p, _3)): p)
            res = dict([(id, (genes, p, ref))
                        for (id, (genes, _, ref)), p in
                        zip(res, obiProb.FDR([p for _, (_, p, _) in res]))])
        return res

    def GetAnnotatedTerms(self, genes, directAnnotationOnly=False,
                          evidenceCodes=None, progressCallback=None):
        """ Return all terms that are annotated by genes with evidenceCodes.
        """
        genes = [genes] if type(genes) == str else genes
        revGenesDict = self.GetGeneNamesTranslator(genes)
        genes = set(revGenesDict.keys())
        evidenceCodes = set(evidenceCodes or evidenceDict.keys())
        annotations = [ann for gene in genes for ann in self.geneAnnotations[gene]
                       if ann.Evidence_Code in evidenceCodes]
        dd = defaultdict(set)
        for ann in annotations:
            dd[ann.GO_ID].add(revGenesDict.get(ann.geneName, ann.geneName))
        if not directAnnotationOnly:
            terms = dd.keys()
            filteredTerms = [term for term in terms if term in self.ontology]
            if len(terms) != len(filteredTerms):
                termDiff = set(terms) - set(filteredTerms)
                warnings.warn(
                    "%s terms in the annotations were not found in the "
                    "ontology." % ",".join(map(repr, termDiff)),
                    UserWarning)

            terms = self.ontology.ExtractSuperGraph(filteredTerms)
            for i, term in enumerate(terms):
                termAnnots = self.GetAllAnnotations(term).intersection(annotations)
##                termAnnots.intersection_update(annotations)
                dd[term].update([revGenesDict.get(ann.geneName, ann.geneName)
                                 for ann in termAnnots])
        return dict(dd)

    def DrawEnrichmentGraph(self, terms, clusterSize, refSize=None,
                            file="graph.png", width=None, height=None,
                            precison=3):
        refSize = len(self.geneNames) if refSize == None else refSize
        sortedterms = sorted(terms.items(), key=lambda term: term[1][1])
        fdr = dict(zip([t[0] for t in sortedterms],
                       obiProb.FDR([t[1][1] for t in sortedterms])))
        termsList = [(term,
                      ((float(len(terms[term][0])) / clusterSize) /
                       (float(terms[term][2]) / refSize)),
                      len(terms[term][0]),
                      terms[term][2],
                      terms[term][1],
                      fdr[term],
                      terms[term][0])
                     for term in terms]

        drawEnrichmentGraph(termsList, file, width, height,
                            ontology=self.ontology, precison=precison)

    def __add__(self, iterable):
        """ Return a new Annotations object with combined annotations
        """
        return Annotations([a for a in self] + [a for a in iterable],
                           ontology=self.ontology)

    def __iadd__(self, iterable):
        """ Add annotations to this instance
        """
        self.extend(iterable)
        return self

    def __contains__(self, item):
        return item in self.annotations

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

    def __getslice__(self, *args):
        return self.annotations.__getslice__(*args)

    def add(self, line):
        """ Add one annotation
        """
        self.AddAnnotation(line)

    def append(self, line):
        """ Add one annotation
        """
        self.AddAnnotation(line)

    def extend(self, lines):
        """ Add multiple annotations
        """
        for line in lines:
            self.AddAnnotation(line)

    def RemapGenes(self, map):
        """
        """
        for gene in map:
            annotations = self.geneAnnotations[gene]
            for ann in annotations:
                for name in map[gene]:
                    ann1 = copy.copy(ann)
                    ann1.geneName = name
                    self.add(ann1)
        self.genematcher = obiGene.GMDirect()
        try:
            del self._geneNames
        except Exception:
            pass
        self.genematcher.set_targets(self.geneNames)

    @staticmethod
    def DownloadAnnotations(org, file, progressCallback=None):
        if isinstance(file, basestring):
            tFile = tarfile.open(file, "w:gz")
        else:
            tFile = file

        tmpDir = os.path.join(orngEnviron.bufferDir, "tmp_go/")
        try:
            os.mkdir(tmpDir)
        except Exception:
            pass
        fileName = "gene_association." + org + ".gz"
        urlretrieve(("http://www.geneontology.org/gene-associations/" +
                     fileName),
                    os.path.join(tmpDir, fileName),
                    progressCallback and __progressCallbackWrapper(progressCallback))
        gzFile = GzipFile(os.path.join(tmpDir, fileName), "r")
        file = open(os.path.join(tmpDir, "gene_association." + org), "w")
        file.writelines(gzFile.readlines())
        file.flush()
        file.close()

        tFile.add(os.path.join(tmpDir, "gene_association." + org),
                  "gene_association")
        annotation = Annotations(os.path.join(tmpDir, "gene_association." + org),
                    genematcher=obiGene.GMDirect(), progressCallback=progressCallback)
        cPickle.dump(annotation.geneNames, open(os.path.join(tmpDir, "gene_names.pickle"), "wb"))
        tFile.add(os.path.join(tmpDir, "gene_names.pickle"), "gene_names.pickle")
        tFile.close()
        os.remove(os.path.join(tmpDir, "gene_association." + org))
        os.remove(os.path.join(tmpDir, "gene_names.pickle"))

    @staticmethod
    def DownloadAnnotationsAtRev(org, rev, filename=None,
                                 progressCallback=None):
        if filename is None:
            filename = os.path.join(default_database_path,
                                    "gene_association.%s@rev%s.tar.gz" %
                                    (code, rev))
        url = ("http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/gene-associations/gene_association.%s.gz?rev=%s" %
               (org, rev))
        url += ";content-type=application%2Fx-gzip"
        r = urllib2.urlopen(url)

        with open(filename + ".part", "wb") as f:
            shutil.copyfileobj(r, f)

        os.rename(filename + ".part", filename)

from .obiTaxonomy import pickled_cache


@pickled_cache(None, [("GO", "taxonomy.pickle"),
                      ("Taxonomy", "ncbi_taxonomy.tar.gz")])
def organism_name_search(name):
    return Annotations.organism_name_search(name)


def filterByPValue(terms, maxPValue=0.1):
    """ Filters the terms by the p-value. Asumes terms is a dict with
    the same structure as returned from GetEnrichedTerms.

    """
    return dict(filter(lambda (k, e): e[1] <= maxPValue, terms.items()))


def filterByFrequency(terms, minF=2):
    """ Filters the terms by the cluster frequency. Asumes terms is
    a dict with the same structure as returned from GetEnrichedTerms.

    """
    return dict(filter(lambda (k, e): len(e[0]) >= minF, terms.items()))


def filterByRefFrequency(terms, minF=4):
    """ Filters the terms by the reference frequency. Asumes terms is
    a dict with the same structure as returned from GetEnrichedTerms.

    """
    return dict(filter(lambda (k, e): e[2] >= minF, terms.items()))


def drawEnrichmentGraph(enriched, file="graph.png", width=None, height=None,
                        header=None, ontology=None, precison=3):
    file = open(file, "wb") if type(file) == str else file
    drawEnrichmentGraph_tostreamMk2(enriched, file, width, height, header,
                                    ontology, precison)


def drawEnrichmentGraph_tostreamMk2(enriched, fh, width, height, header=None,
                                    ontology=None, precison=4):
    ontology = ontology if ontology else Ontology()
    header = header if header else ["List", "Total", "p-value", "FDR",
                                    "Names", "Genes"]
    GOTerms = dict([(t[0], t) for t in enriched if t[0] in ontology])

    def getParents(term):
        parents = ontology.ExtractSuperGraph([term])
        parents = [id for id in parents if id in GOTerms and id != term]
        c = reduce(set.union, [set(ontology.ExtractSuperGraph([id])) -
                               set([id]) for id in parents], set())
        parents = [t for t in parents if t not in c]
        return parents
    parents = dict([(term, getParents(term)) for term in GOTerms])
    # print "Parentes", parents

    def getChildren(term):
        return [id for id in GOTerms if term in parents[id]]
    topLevelTerms = [id for id in parents if not parents[id]]
    # print "Top level terms", topLevelTerms
    termsList = []
    fmt = "%" + ".%if" % precison

    def collect(term, parent):
        termsList.append(GOTerms[term][1:4] + \
                         (fmt % GOTerms[term][4],
                          fmt % GOTerms[term][5],
                          ontology[term].name,
                          ", ".join(GOTerms[term][6])) + (parent,))
        parent = len(termsList) - 1
        for c in getChildren(term):
            collect(c, parent)

    for topTerm in topLevelTerms:
        collect(topTerm, None)
    for entry in enriched:
        if entry[0] not in ontology:
            termsList.append(entry[1:4] + \
                             (fmt % entry[4],
                              fmt % entry[5],
                              entry[0],
                              ", ".join(entry[6])) + (None,))

    drawEnrichmentGraphPIL_tostream(termsList, header, fh, width, height)


def drawEnrichmentGraphPIL_tostream(termsList, headers, fh, width=None, height=None):
    from PIL import Image, ImageDraw, ImageFont
    backgroundColor = (255, 255, 255)
    textColor = (0, 0, 0)
    graphColor = (0, 0, 255)
    fontSize = height == None and 12 or (height - 60) / len(termsList)
    font = ImageFont.load_default()
    try:
        font = ImageFont.truetype("arial.ttf", fontSize)
    except:
        pass
    getMaxTextHeightHint = lambda l: max([font.getsize(t)[1] for t in l])
    getMaxTextWidthHint = lambda l: max([font.getsize(t)[0] for t in l])
    maxFoldWidth = width != None and min(150, width / 6) or 150
    maxFoldEnrichment = max([t[0] for t in termsList])
    foldNormalizationFactor = float(maxFoldWidth) / maxFoldEnrichment
    foldWidths = [int(foldNormalizationFactor * term[0]) for term in termsList]
    treeStep = 10
    treeWidth = {}
    for i, term in enumerate(termsList):
        treeWidth[i] = (term[-1] == None and 1 or treeWidth[term[-1]] + 1)
    treeStep = width != None and min(treeStep, width / (6 * max(treeWidth.values())) or 2) or treeStep
    treeWidth = [w * treeStep + foldWidths[i] for i, w in treeWidth.items()]
    treeWidth = max(treeWidth) - maxFoldWidth
    verticalMargin = 10
    horizontalMargin = 10
##    print verticalMargin, maxFoldWidth, treeWidth
##    treeWidth = 100
    firstColumnStart = verticalMargin + maxFoldWidth + treeWidth + 10
    secondColumnStart = firstColumnStart + getMaxTextWidthHint([str(t[1]) for t in termsList] + [headers[0]]) + 2
    thirdColumnStart = secondColumnStart + getMaxTextWidthHint([str(t[2]) for t in termsList] + [headers[1]]) + 2
    fourthColumnStart = thirdColumnStart + getMaxTextWidthHint([str(t[3]) for t in termsList] + [headers[2]]) + 2
    fifthColumnStart = fourthColumnStart + getMaxTextWidthHint([str(t[4]) for t in termsList] + [headers[3]]) + 4
##    maxAnnotationTextWidth = width==None and getMaxTextWidthHint([str(t[4]) for t in termsList]+["Annotation"]) or (width - fourthColumnStart - verticalMargin) * 2 / 3
    maxAnnotationTextWidth = width == None and getMaxTextWidthHint([str(t[5]) for t in termsList] + [headers[4]]) or max((width - fifthColumnStart - verticalMargin) * 2 / 3, getMaxTextWidthHint([t[5] for t in termsList] + [headers[4]]))
    sixthColumnStart = fifthColumnStart + maxAnnotationTextWidth + 4
    maxGenesTextWidth = width == None and getMaxTextWidthHint([str(t[6]) for t in termsList] + [headers[5]]) or (width - fifthColumnStart - verticalMargin) / 3

    legendHeight = font.getsize("1234567890")[1] * 2
    termHeight = font.getsize("A")[1]
##    print fourthColumnStart, maxAnnotationTextWidth, verticalMargin
    width = sixthColumnStart + maxGenesTextWidth + verticalMargin
    height = len(termsList) * termHeight + 2 * (legendHeight + horizontalMargin)

    image = Image.new("RGB", (width, height), backgroundColor)
    draw = ImageDraw.Draw(image)

    def truncText(text, maxWidth, append=""):
        # print getMaxTextWidthHint([text]), maxAnnotationTextWidth
        if getMaxTextWidthHint([text]) > maxWidth:
            while getMaxTextWidthHint([text + "..." + append]) > maxWidth and text:
                text = text[:-1]
            if text:
                text = text + "..." + append
            else:
                text = append
        return text
    currentY = horizontalMargin + legendHeight
    connectAtX = {}
    for i, term in enumerate(termsList):
        draw.line([(verticalMargin, currentY + termHeight / 2), (verticalMargin + foldWidths[i], currentY + termHeight / 2)], width=termHeight - 2, fill=graphColor)
        draw.text((firstColumnStart, currentY), str(term[1]), font=font, fill=textColor)
        draw.text((secondColumnStart, currentY), str(term[2]), font=font, fill=textColor)
        draw.text((thirdColumnStart, currentY), str(term[3]), font=font, fill=textColor)
        draw.text((fourthColumnStart, currentY), str(term[4]), font=font, fill=textColor)
##        annotText = width!=None and truncText(str(term[5]), maxAnnotationTextWidth, str(term[5])) or str(term[4])
        annotText = width != None and truncText(str(term[5]), maxAnnotationTextWidth)
        draw.text((fifthColumnStart, currentY), annotText, font=font, fill=textColor)
        genesText = width != None and truncText(str(term[6]), maxGenesTextWidth) or str(term[6])
        draw.text((sixthColumnStart, currentY), genesText, font=font, fill=textColor)
        lineEnd = term[-1] == None and firstColumnStart - 10 or connectAtX[term[-1]]
        draw.line([(verticalMargin + foldWidths[i] + 1, currentY + termHeight / 2), (lineEnd, currentY + termHeight / 2)], width=1, fill=textColor)
        if term[-1] != None:
            draw.line([(lineEnd, currentY + termHeight / 2), (lineEnd, currentY + termHeight / 2 - termHeight * (i - term[-1]))], width=1, fill=textColor)
        connectAtX[i] = lineEnd - treeStep
        currentY += termHeight

    currentY = horizontalMargin
    draw.text((firstColumnStart, currentY), headers[0], font=font, fill=textColor)
    draw.text((secondColumnStart, currentY), headers[1], font=font, fill=textColor)
    draw.text((thirdColumnStart, currentY), headers[2], font=font, fill=textColor)
    draw.text((fourthColumnStart, currentY), headers[3], font=font, fill=textColor)
    draw.text((fifthColumnStart, currentY), headers[4], font=font, fill=textColor)
    draw.text((sixthColumnStart, currentY), headers[5], font=font, fill=textColor)

    horizontalMargin = 0
    # draw.line([(verticalMargin, height - horizontalMargin - legendHeight), (verticalMargin + maxFoldWidth, height - horizontalMargin - legendHeight)], width=1, fill=textColor)
    draw.line([(verticalMargin, horizontalMargin + legendHeight), (verticalMargin + maxFoldWidth, horizontalMargin + legendHeight)], width=1, fill=textColor)
    maxLabelWidth = getMaxTextWidthHint([" " + str(i) for i in range(int(maxFoldEnrichment + 1))])
    numOfLegendLabels = max(int(maxFoldWidth / maxLabelWidth), 2)
    for i in range(numOfLegendLabels + 1):
        # draw.line([(verticalMargin + i*maxFoldWidth/10, height - horizontalMargin - legendHeight/2), (verticalMargin + i*maxFoldWidth/10, height - horizontalMargin - legendHeight)], width=1, fill=textColor)
        # draw.text((verticalMargin + i*maxFoldWidth/10 - font.getsize(str(i))[0]/2, height - horizontalMargin - legendHeight/2), str(i), font=font, fill=textColor)

        label = str(int(i * maxFoldEnrichment / numOfLegendLabels))
        draw.line([(verticalMargin + i * maxFoldWidth / numOfLegendLabels, horizontalMargin + legendHeight / 2), (verticalMargin + i * maxFoldWidth / numOfLegendLabels, horizontalMargin + legendHeight)], width=1, fill=textColor)
        draw.text((verticalMargin + i * maxFoldWidth / numOfLegendLabels - font.getsize(label)[0] / 2, horizontalMargin), label, font=font, fill=textColor)

    image.save(fh)


def drawEnrichmentGraphPylab_tostream(termsList, headers, fh, width=None, height=None, show=True):
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle

    maxFoldWidth = width != None and min(150, width / 6) or 150
    maxFoldEnrichment = max([t[0] for t in termsList])
    foldNormalizationFactor = float(maxFoldWidth) / maxFoldEnrichment
##    foldWidths = [int(foldNormalizationFactor*term[0]) for term in termsList]
    foldWidths = [term[0] for term in termsList]
    treeStep = maxFoldEnrichment * 0.05
    treeWidth = {}

    for i, term in enumerate(termsList):
        treeWidth[i] = (term[-1] == None and treeStep or treeWidth[term[-1]] + treeStep)
    maxTreeWidth = max(treeWidth)

    connectAt = {}
    cellText = []
    axes1 = plt.axes([0.1, 0.1, 0.2, 0.8])
    for i, line in enumerate(termsList):
        enrichment, n, m, p_val, fdr_val, name, genes, parent = line
        r = Rectangle((0, len(termsList) - i - 0.4), enrichment, 0.8)
        plt.gca().add_patch(r)
        plt.plot([enrichment, connectAt.get(parent, maxFoldEnrichment + maxTreeWidth)], [len(termsList) - i, len(termsList) - i], color="black")
        connectAt[i] = connectAt.get(parent, maxFoldEnrichment + maxTreeWidth) - treeStep
        if parent != None:
            plt.plot([connectAt.get(parent)] * 2, [len(termsList) - i, len(termsList) - parent], color="black")
        cellText.append((str(n), str(m), p_val, fdr_val, name, genes))

##    from Orange.orng.orngClustering import TableTextLayout
##    text = TableTextLayout((maxFoldEnrichment*1.1, len(termsList)), cellText)
    from Orange.orng.orngClustering import TablePlot
    if True:
        axes2 = plt.axes([0.3, 0.1, 0.6, 0.8], sharey=axes1)
        axes2.set_axis_off()
        table = TablePlot((0, len(termsList)), axes=plt.gca())
        for i, line in enumerate(cellText):
            for j, text in enumerate(line):
                table.add_cell(i, j, width=len(text), height=1, text=text, loc="left", edgecolor="w", facecolor="w")

        table.set_figure(plt.gcf())
        plt.gca().add_artist(table)
        plt.gca()._set_artist_props(table)
##    plt.text(3, 3, "\n".join(["\t".join(text) for text in cellText]))

##    table = plt.table(cellText=cellText, colLabels=headers, loc="right")
##    table.set_transform(plt.gca().transData)
##
##    table.set_xy(20,20)
    plt.show()


class Taxonomy(object):
    """Maps NCBI taxonomy ids to corresponding GO organism codes
    """
    common_org_map = {"297284": "9913", "30523": "9913",  # Bos taurus
                      "5782": "352472", "44689": "352472", "366501": "352472",  # Dictyostelium discoideum
                      "83333": "562",  # Escherichia coli
                      "52545": "4530", "4532": "4530", "65489": "4530", "4533": "4530", "77588": "4530", "29689": "4530",
                      "4538": "4530", "40148": "4530", "29690": "4530", "110450": "4530", "4534": "4530", "83309": "4530",
                      "4528": "4530", "127571": "4530", "40149": "4530", "83307": "4530", "63629": "4530", "4536": "4530",
                      "4535": "4530", "4537": "4530", "65491": "4530", "83308": "4530", "4529": "4530", "4530": "4530",
                      "39946": "4530", "39947": "4530", "110451": "4530", "364100": "4530", "364099": "4530", "4539": "4530",
                      }
    code_map = {"3702": "tair",  # Arabidopsis thaliana
                "9913": "goa_cow",  # Bos taurus
                "6239": "wb",  # Caenorhabditis elegans
                "3055": None,  # Chlamydomonas reinhardtii
                "7955": "zfin",  # Danio rerio (zebrafish)
                "352472": "dictyBase",  # Dictyostelium discoideum
                "7227": "fb",  # Drosophila melanogaster
                "562": "ecocyc",  # Escherichia coli
                "11103": None,  # Hepatitis C virus
                "9606": "goa_human",  # Homo sapiens
                "10090": "mgi",  # Mus musculus
                "2104": None,  # Mycoplasma pneumoniae
                "4530": "gramene_oryza",  # Oryza sativa
                "5833": "GeneDB_Pfalciparum",  # Plasmodium falciparum
                "4754": None,  # Pneumocystis carinii
                "10116": "rgd",  # Rattus norvegicus
                "4932": "sgd",  # Saccharomyces cerevisiae
                "4896": "GeneDB_Spombe",  # Schizosaccharomyces pombe
                "31033": None,  # Takifugu rubripes
                "8355": None,  # Xenopus laevis
                "4577": None  # Zea mays
                }
    version = 1
    __shared_state = {"tax": None}

    def __init__(self):
        self.__dict__ = self.__shared_state
        if not self.tax:
            path = orngServerFiles.localpath_download("GO", "taxonomy.pickle")
            if os.path.isfile(path):
                self.tax = cPickle.load(open(path, "rb"))
            else:
                orngServerFiles.download("GO", "taxonomy.pickle")
                self.tax = cPickle.load(open(path, "rb"))

    def __getitem__(self, key):
        key = self.common_org_map.get(key, key)
        return self.code_map[key]

    def keys(self):
        return list(set(self.common_org_map.keys() + self.code_map.keys()))


def from_taxid(id):
    """ Return a set of GO organism codes that correspond to NCBI taxonomy id
    """
    return Taxonomy()[id]


def to_taxid(db_code):
    """ Return a set of NCBI taxonomy ids from db_code GO organism annotations
    """
    r = [key for key, val in Taxonomy().code_map.items() if db_code == val]
    return set(r)


class __progressCallbackWrapper:
    def __init__(self, callback):
        self.callback = callback

    def __call__(self, bCount, bSize, fSize):
        fSize = 10000000 if fSize == -1 else fSize
        self.callback(100 * bCount * bSize / fSize)

from .obiGenomicsUpdate import Update as UpdateBase


class Update(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path or getDataDir(), progressCallback)

    def CheckModified(self, addr, date=None):
        return date < self.GetLastModified(addr) if date else True

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
                for taxId in [t.split(":")[-1] for t in taxons if "|" not in t]:  # exclude taxons with cardinality 2
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

    a.GetEnrichedTerms(sorted(a.geneNames)[:100])  # , progressCallback=_print)
##    profile.runctx("a.GetEnrichedTerms(sorted(a.geneNames)[:100])", {"a":a}, {})
    a.GetEnrichedTerms(sorted(a.geneNames)[:100])  # , progressCallback=_print)
    d1 = a.GetEnrichedTerms(sorted(a.geneNames)[:1000])  # , progressCallback=_print)

##    print a.GetEnrichedTerms(sorted(a.geneNames)[:100])#, progressCallback=_print)


def _test2():
    o = Ontology()
    a = Annotations("human", ontology=o)
    clusterGenes = sorted(a.geneNames)[:100]
    for i in range(10):
        terms = a.GetEnrichedTerms(sorted(a.geneNames)[:200], aspect=["P"])
        a.GetEnrichedTerms(sorted(a.geneNames)[:200], aspect=["C"])
        a.GetEnrichedTerms(sorted(a.geneNames)[:200], aspect=["F"])
        print i
#    a.DrawEnrichmentGraph(filterByPValue(terms), len(clusterGenes), len(a.geneNames))

##    drawEnrichmentGraph([("bal", 1.0, 5, 6, 0.1, 0.4, ["vv"]),
##                        ("GO:0019079", 0.5, 5, 6, 0.1, 0.4, ["cc", "bb"]),
##                        ("GO:0022415", 0.4, 5, 7, 0.11, 0.4, ["cc1", "bb"])], open("graph.png", "wb"), None, None)


def _test3():
    o = Ontology()
    a = Annotations("sgd", ontology=o)
##    a = Annotations(list(a)[3:len(a)/3], ontology=o)
    clusterGenes = sorted(a.geneNames)[:1] + sorted(a.geneNames)[-1:]
##    clusterGenes = [g + "_" + str(i%5) for g in sorted(a.geneNames)[:2]]
    exonMap = dict([(gene, [gene + "_E%i" % i for i in range(10)]) for gene in a.geneNames])
    a.RemapGenes(exonMap)
##    o.reverseAliasMapper = o.aliasMapper = {}
    terms = a.GetEnrichedTerms(exonMap.values()[0][:2] + exonMap.values()[-1][2:])
##    terms = a.GetEnrichedTerms(clusterGenes)
    print terms
##    a.DrawEnrichmentGraph(filterByPValue(terms), len(clusterGenes), len(a.geneNames))
    a.DrawEnrichmentGraph(filterByPValue(terms, maxPValue=0.1), len(clusterGenes), len(a.geneNames))

if __name__ == "__main__":
    _test2()
