"""
Gene Ontology (:mod:`go`)
=========================

"""

from __future__ import absolute_import

import os
import tarfile
import gzip
import re
import cPickle
import shutil
import urllib2
import warnings

from gzip import GzipFile
from collections import defaultdict
from operator import attrgetter

from Orange.utils import \
    deprecated_keywords, deprecated_members, progress_bar_milestones, \
    serverfiles, environ

from .utils import stats

from . import gene as obiGene, taxonomy as obiTaxonomy

default_database_path = os.path.join(serverfiles.localpath(), "GO")


_CVS_REVISION_RE = re.compile(r"^(rev)?(\d+\.\d+)+$")

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
    "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By",
    # GAF v2.0
    "Annotation_Extension", "Gene_Product_Form_ID"
]

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


@deprecated_members({"ParseStanza": "parse_stanza",
                     "GetRelatedObjcts": "related_objects",
                     "relatedTo": "related_to"})
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
        self.related_to = set()
        if stanza:
            self.parse_stanza(stanza)

    def parse_stanza(self, stanza):
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
        self.related = set(self.related_objects())
        self.__dict__.update(self.values)
        if "def" in self.__dict__:
            self.__dict__["def_"] = self.__dict__["def"]

    def related_objects(self):
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
        for typeId, id in self.related_to:
            yield (typeId, self.ontology[id])


class Term(OBOObject):
    pass


class Typedef(OBOObject):
    pass


class Instance(OBOObject):
    pass


@deprecated_members(
    {"ParseFile": "parse_file", "slimsSubset": "slims_subset",
     "GetDefinedSlimsSubsets": "defined_slim_subsets",
     "SetSlimSubsets": "set_slim_subsets",
     "GetSlimsSubset": "named_slims_subset",
     "GetSlimTerms": "slims_for_term",
     "ExtractSuperGraph": "extract_super_graph",
     "ExtractSubGraph": "extract_sub_graph",
     "GetTermDepth": "term_depth",
     "aliasMapper": "alias_mapper",
     "reverseAliasMapper": "reverse_alias_mapper"},
    wrap_methods=[])
class Ontology(object):
    """
    :class:`Ontology` is the class representing a gene ontology.

    :param str filename:
        A filename of an .obo formated file.
    :param progress_callback:
        Optional `float -> None` function.
    :param str rev:
        An CVS revision specifier (see `GO web CVS interface
        <http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/ontology/>`_)

    Example:

        >>> # Load the current ontology (downloading it if necessary)
        >>> ontology = Ontology()
        >>> # Load the ontology at the specified CVS revision.
        >>> ontology = Ontology(rev="5.2092")

    """
    version = 1

    @deprecated_keywords({"progressCallback": "progress_callback"})
    def __init__(self, filename=None, progress_callback=None, rev=None):
        self.terms = {}
        self.typedefs = {}
        self.instances = {}
        self.slims_subset = set()
        self.alias_mapper = {}
        self.reverse_alias_mapper = defaultdict(set)
        self.header = ""

        if filename is not None:
            self.parse_file(filename, progress_callback)
        elif rev is not None:
            if not _CVS_REVISION_RE.match(rev):
                raise ValueError("Invalid revision format.")
            if rev.startswith("rev"):
                rev = rev[3:]
            pc = lambda v: progress_callback(v / 2.0) \
                 if progress_callback else None
            filename = os.path.join(default_database_path,
                                    "gene_ontology_edit@rev%s.obo" % rev)
            if not os.path.exists(filename):
                self.download_ontology_at_rev(rev, filename, pc)
            self.parse_file(filename,
                            lambda v: progress_callback(v / 2.0 + 50)
                            if progress_callback else None)
        else:
            filename = serverfiles.localpath_download(
                "GO", "gene_ontology_edit.obo.tar.gz"
            )
            self.parse_file(filename, progress_callback)

    @classmethod
    @deprecated_keywords({"progressCallback": "progress_callback"})
    def load(cls, progress_callback=None):
        """ A class method that tries to load the ontology file from
        default_database_path. It looks for a filename starting with
        'gene_ontology'. If not found it will download it.

        """
        filename = os.path.join(default_database_path,
                                "gene_ontology_edit.obo.tar.gz")
        if not os.path.isfile(filename) and not os.path.isdir(filename):
            serverfiles.download("GO", "gene_ontology_edit.obo.tar.gz")

        return cls(filename, progress_callback=progress_callback)

    Load = load

    @deprecated_keywords({"progressCallback": "progress_callback"})
    def parse_file(self, file, progress_callback=None):
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

        milestones = progress_bar_milestones(len(data), 90)
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
            if progress_callback and i in milestones:
                progress_callback(90.0 * i / len(data))

        self.alias_mapper = {}
        self.reverse_alias_mapper = defaultdict(set)
        milestones = progress_bar_milestones(len(self.terms), 10)
        for i, (id, term) in enumerate(self.terms.iteritems()):
            for typeId, parent in term.related:
                self.terms[parent].related_to.add((typeId, id))
            try:
                self.alias_mapper.update([(alt_id, id)
                                          for alt_id in term.alt_id])
                self.reverse_alias_mapper[id].union_update(term.alt_id)
            except AttributeError:
                pass
            if progress_callback and i in milestones:
                progress_callback(90.0 + 10.0 * i / len(self.terms))

    def defined_slims_subsets(self):
        """
        Return a list of defined subsets in the ontology.

        :rtype: list-of-str

        """
        return [line.split()[1] for line in self.header.splitlines()
                if line.startswith("subsetdef:")]

    def named_slims_subset(self, subset):
        """
        Return all term IDs in a named `subset`.

        :param str subset: A string naming a subset in the ontology.
        :rtype: list-of-str

        .. seealso:: :func:`defined_slims_subsets`

        """
        return [id for id, term in self.terms.items()
                if subset in getattr(term, "subset", set())]

    def set_slims_subset(self, subset):
        """
        Set the `slims_subset` term subset to `subset`.

        :param set subset: A subset of GO term IDs.

        `subset` may also be a string, in which case the call is equivalent
        to ``ont.set_slims_subsets(ont.named_slims_subset(subset))``

        """
        if isinstance(subset, basestring):
            self.slims_subset = set(self.named_slims_subset(subset))
        else:
            self.slims_subset = set(subset)

    def slims_for_term(self, term):
        """
        Return a list of slim term IDs for `term`.

        This is a list of `most specific` slim terms to which `term` belongs.

        :param str term: Term ID.

        """
        queue = set([term])
        visited = set()
        slims = set()
        while queue:
            term = queue.pop()
            visited.add(term)
            if term in self.slims_subset:
                slims.add(term)
            else:
                queue.update(set(tid for _, tid in self[term].related) -
                             visited)
        return slims

    def extract_super_graph(self, terms):
        """
        Return all super terms of `terms` up to the most general one.

        :param list terms: A list of term IDs.

        """
        terms = [terms] if isinstance(terms, basestring) else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(set(tid for _, tid in self[term].related) -
                         visited)
        return visited

    def extract_sub_graph(self, terms):
        """
        Return all sub terms of `terms`.

        :param list terms: A list of term IDs.

        """
        terms = [terms] if type(terms) == str else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(set(tid for _, tid in self[term].related_to) -
                         visited)
        return visited

    def term_depth(self, term, cache_={}):
        """
        Return the minimum depth of a `term`.

        (length of the shortest path to this term from the top level term).

        """
        if term not in cache_:
            cache_[term] = min([self.term_depth(parent) + 1
                                for _, parent in self[term].related] or
                               [1])
        return cache_[term]

    def __getitem__(self, termid):
        """
        Return a :class:`Term` object with `termid`.

        :param str term: An id of a 'Term' in the ontology.
        :rtype: :class:`Term`

        """
        if termid in self.terms:
            return self.terms[termid]
        elif termid in self.alias_mapper:
            return self.terms[self.alias_mapper[termid]]
        else:
            raise KeyError(termid)

    def __iter__(self):
        """
        Iterate over all term ids in ontology.
        """
        return iter(self.terms)

    def __len__(self):
        """
        Return number of terms in ontology.
        """
        return len(self.terms)

    def __contains__(self, termid):
        """
        Return `True` if a term with `termid` is present in the ontology.
        """
        return termid in self.terms or termid in self.alias_mapper

    @staticmethod
    @deprecated_keywords({"progressCallback": "progress_callback"})
    def download_ontology(file, progress_callback=None):
        tFile = tarfile.open(file, "w:gz") if isinstance(file, basestring) \
                else file
        tmpDir = os.path.join(environ.buffer_dir, "tmp_go/")
        try:
            os.mkdir(tmpDir)
        except Exception:
            pass

        from Orange.utils import wget
        wget("http://www.geneontology.org/ontology/gene_ontology_edit.obo",
             tmpDir,
             progress=progress_callback)

        tFile.add(os.path.join(tmpDir, "gene_ontology_edit.obo"),
                  "gene_ontology_edit.obo")
        tFile.close()
        os.remove(os.path.join(tmpDir, "gene_ontology_edit.obo"))

    DownloadOntology = download_ontology

    @staticmethod
    @deprecated_keywords({"progressCallback": "progress_callback"})
    def download_ontology_at_rev(rev, filename=None, progress_callback=None):
        url = "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/ontology/gene_ontology_edit.obo?rev=%s" % rev
        url += ";content-type=text%2Fplain"
        if filename is None:
            filename = os.path.join(default_database_path,
                                    "gene_ontology_edit@rev%s.obo" % rev)
        r = urllib2.urlopen(url)

        with open(filename + ".part", "wb") as f:
            shutil.copyfileobj(r, f)

        os.rename(filename + ".part", filename)

    DownloadOntologyAtRev = download_ontology_at_rev


from collections import namedtuple

_AnnotationRecordBase = namedtuple(
    "AnnotationRecord",
    annotationFields
)


class AnnotationRecord(_AnnotationRecordBase):
    """
    An annotation record mapping a gene to a term.

    See http://geneontology.org/GO.format.gaf-2_0.shtml for description
    if individual fields.

    """
    def __new__(cls, *args):
        if len(args) == 1 and isinstance(args[0], basestring):
            args = map(intern, args[0].split("\t"))
        return super(AnnotationRecord, cls).__new__(cls, *args)

    @classmethod
    def from_string(cls, string):
        """
        Create an instance from a line in a annotations (GAF 2.0 format) file.
        """
        return AnnotationRecord._make(map(intern, string.split("\t")))

    gene_name = property(
        attrgetter("DB_Object_Symbol"),
        doc="Alias for DB_Object_Symbol"
    )
    geneName = gene_name
    GOId = property(
        attrgetter("GO_ID"),
        doc="Alias for GO_ID"
    )
    go_id = GOId

    evidence = property(
        attrgetter("Evidence_Code"),
        doc="Alias for Evidence_Code"
    )
    aspect = property(
        attrgetter("Aspect"),
        doc="Alias for Aspect"
    )

    @property
    def alias(self):
        return list(map(intern, self.DB_Object_Synonym.split("|")))


@deprecated_members(
    {"GetOntology": "get_ontology", "SetOntology": "set_ontology",
     "ParseFile": "parse_file", "AddAnnotation": "add_annotation",
     "GetGeneNamesTranslator": "get_gene_names_translator",
     "GetAllAnnotations": "get_all_annotations",
     "GetAllGenes": "get_all_genes",
     "GetEnrichedTerms": "get_enriched_terms",
     "GetAnnotatedTerms": "get_annotated_terms",
     "DrawEnrichmentGraph": "draw_enrichment_graph",
     "geneNamesDict": "gene_names_dict", "geneNames": "gene_names",
     "aliasMapper": "alias_mapper",
     "allAnnotations": "all_annotations",
     "geneAnnotations": "gene_annotations",
     "termAnnotations": "term_annotations"},
    wrap_methods=[])
class Annotations(object):
    """
    :class:`Annotations` object holds the annotations.

    :param str filename_or_org:
        A filename of a GAF formated annotations file (e.g.
        gene_annotations.goa_human) or an organism specifier (e.g.
        ``'goa_human'`` or ``'9606'``). In the later case the annotations
        for that organism will be loaded.

    :param ontology: :class:`Ontology` object for annotations
    :type ontology: :class:`Ontology`

    :param str rev:
        An optional CVS revision string. If the `filename_or_org` is given an
        organism code the annotations will be retrieved for that revision
        (see `GO web CVS
        <http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/>`_)

    """
    version = 2

    @deprecated_keywords({"progressCallback": "progress_callback"})
    def __init__(self, filename_or_organism=None, ontology=None, genematcher=None,
                 progress_callback=None, rev=None):
        self.ontology = ontology

        #: A dictionary mapping a gene name (DB_Object_Symbol) to a
        #: set of all annotations of that gene.
        self.gene_annotations = defaultdict(list)

        #: A dictionary mapping a GO term id to a set of annotations that
        #: are directly annotated to that term
        self.term_anotations = defaultdict(list)

        self.all_annotations = defaultdict(list)

        self._gene_names = None
        self._gene_names_dict = None
        self._alias_mapper = None

        #: A list of all :class:`AnnotationRecords` instances.
        self.annotations = []
        self.header = ""
        self.genematcher = genematcher
        self.taxid = None

        if type(filename_or_organism) in [list, set, dict, Annotations]:
            for ann in filename_or_organism:
                self.add_annotation(ann)
            if type(filename_or_organism, Annotations):
                self.taxid = filename_or_organism.taxid

        elif isinstance(filename_or_organism, basestring) and \
                os.path.exists(filename_or_organism):
            self.parse_file(filename_or_organism, progress_callback)

        elif isinstance(filename_or_organism, basestring):
            # Assuming organism code/name
            if rev is not None:
                if not _CVS_REVISION_RE.match(rev):
                    raise ValueError("Invalid revision format")

                if rev.startswith("rev"):
                    rev = rev[3:]
                code = self.organism_name_search(filename_or_organism)
                filename = os.path.join(default_database_path,
                                        "gene_association.%s@rev%s.tar.gz" %
                                        (code, rev))

                if not os.path.exists(filename):
                    self.DownloadAnnotationsAtRev(
                        code, rev, filename, progress_callback)

                self.parse_file(filename, progress_callback)
                self.taxid = to_taxid(code).pop()
            else:
                a = self.Load(filename_or_organism, ontology, genematcher, progress_callback)
                self.__dict__ = a.__dict__
                self.taxid = to_taxid(organism_name_search(filename_or_organism)).pop()
        elif filename_or_organism is not None:
            self.parse_file(filename_or_organism, progress_callback)

        if not self.genematcher and self.taxid:
            matchers = [obiGene.GMGO(self.taxid)]
            if self.taxid == "352472":
                matchers.extend(
                    [obiGene.GMDicty(),
                     [obiGene.GMGO(self.taxid), obiGene.GMDicty()]]
                )

            self.genematcher = obiGene.matcher(matchers)

        if self.genematcher:
            self.genematcher.set_targets(self.gene_names)

    @classmethod
    def organism_name_search(cls, org):
        ids = to_taxid(org)
        from . import taxonomy as tax
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
        serverfiles.localpath_download(
            "GO", "gene_association.%s.tar.gz" % name)
        return ("v%i." % cls.version) + serverfiles.info("GO",
                        "gene_association.%s.tar.gz" % name)["datetime"]

    def set_ontology(self, ontology):
        """Set the ontology to use in the annotations mapping.
        """
        self.all_annotations = defaultdict(list)
        self._ontology = ontology

    def get_ontology(self):
        return self._ontology

    ontology = property(get_ontology, set_ontology,
                        doc=":class:`Ontology` object for annotations.")

    def _ensure_ontology(self):
        if self.ontology is None:
            self.ontology = Ontology()

    @classmethod
    @deprecated_keywords({"progressCallback": "progress_callback"})
    def load(cls, org, ontology=None, genematcher=None,
             progress_callback=None):
        """A class method that tries to load the association file for the
        given organism from default_database_path.
        """
        code = organism_name_search(org)

        filename = "gene_association.%s.tar.gz" % code

        path = serverfiles.localpath("GO", filename)

        if not os.path.exists(path):
            sf = serverfiles.ServerFiles()
            available = sf.listfiles("GO")
            if filename not in available:
                raise obiTaxonomy.UnknownSpeciesIdentifier(org + str(code))
            serverfiles.download("GO", filename)

        return cls(path, ontology=ontology, genematcher=genematcher,
                   progress_callback=progress_callback)

    Load = load

    @deprecated_keywords({"progressCallback": "progress_callback"})
    def parse_file(self, file, progress_callback=None):
        """Parse and load the annotations from file.

        `file` can be:
            - a tarball containing the association file named gene_association
            - a directory name containing the association file named
              gene_association
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

        milestones = progress_bar_milestones(len(lines), 100)
        for i, line in enumerate(lines):
            if line.startswith("!"):
                self.header = self.header + line + "\n"
                continue

            a = AnnotationRecord.from_string(line)
            self.add_annotation(a)

            if progress_callback and i in milestones:
                progress_callback(100.0 * i / len(lines))

    def add_annotation(self, a):
        """Add a single :class:`AnotationRecord` instance to this object.
        """
        if not isinstance(a, AnnotationRecord):
            a = AnnotationRecord(a)
        if not a.geneName or not a.GOId or a.Qualifier == "NOT":
            return

        self.gene_annotations[a.geneName].append(a)
        self.annotations.append(a)
        self.term_anotations[a.GOId].append(a)
        self.all_annotations = defaultdict(list)

        self._gene_names_dict = None
        self._gene_names = None
        self._alias_mapper = None

    @property
    def gene_names_dict(self):
        if self._gene_names_dict is None:
            self._gene_names_dict = defaultdict(set)
            for alias, name in self.alias_mapper.iteritems():
                self._gene_names_dict[name].add(alias)
        return self._gene_names_dict

    @property
    def gene_names(self):
        if self._gene_names is None:
            self._gene_names = set([ann.geneName for ann in self.annotations])
        return self._gene_names

    @property
    def alias_mapper(self):
        if self._alias_mapper is None:
            self._alias_mapper = {}
            for ann in self.annotations:
                self._alias_mapper.update([(alias, ann.geneName)
                                           for alias in ann.alias +
                                            [ann.geneName, ann.DB_Object_ID]])
        return self._alias_mapper

    def get_gene_names_translator(self, genes):
        """ Return a dictionary mapping canonical names (DB_Object_Symbol)
        to `genes`.

        """
        def alias(gene):
            if self.genematcher:
                return self.genematcher.umatch(gene)
            else:
                return (gene if gene in self.gene_names
                        else self.alias_mapper.get(gene, None))

        return dict([(alias(gene), gene) for gene in genes if alias(gene)])

    def _collect_annotations(self, id, visited):
        """ Recursive function collects and caches all annotations for id
        """
        if id not in self.all_annotations and id not in visited:
            if id in self.ontology.reverse_alias_mapper:
                annotations = [self.term_anotations.get(alt_id, [])
                               for alt_id in
                               self.ontology.reverse_alias_mapper[id]] + \
                              [self.term_anotations[id]]
            else:
                ## annotations for this term alone
                annotations = [self.term_anotations[id]]
            visited.add(id)
            for typeId, child in self.ontology[id].related_to:
                aa = self._collect_annotations(child, visited)
                if type(aa) == set:
                    ## if it was already reduced in get_all_annotations
                    annotations.append(aa)
                else:
                    annotations.extend(aa)
            self.all_annotations[id] = annotations
        return self.all_annotations[id]

    _CollectAnnotations = _collect_annotations

    def get_all_annotations(self, id):
        """ Return a set of all annotations (instances of :obj:`AnnotationRecord`)
        for GO term `id` and all it's subterms.

        :param str id: GO term id

        """
        self._ensure_ontology()
        id = self.ontology.alias_mapper.get(id, id)
        if id not in self.all_annotations or \
                type(self.all_annotations[id]) == list:
            annot_set = set()
            for annots in self._collect_annotations(id, set()):
                annot_set.update(annots)
            self.all_annotations[id] = annot_set
        return self.all_annotations[id]

    @deprecated_keywords({"evidenceCodes": "evidence_codes"})
    def get_all_genes(self, id, evidence_codes=None):
        """ Return a list of genes annotated by specified `evidence_codes`
        to GO term 'id' and all it's subterms."

        :param str id: GO term id

        :param list-of-strings evidence_codes:
            List of evidence codes to consider when matching annotations
            to terms.

        """
        evidence_codes = set(evidence_codes or evidenceDict.keys())
        annotations = self.get_all_annotations(id)
        return list(set([ann.geneName for ann in annotations
                         if ann.Evidence_Code in evidence_codes]))

    @deprecated_keywords({
        "evidenceCodes": "evidence_codes", "slimsOnly": "slims_only",
        "useFDR": "use_fdr", "progressCallback": "progress_callback"})
    def get_enriched_terms(self, genes, reference=None, evidence_codes=None,
                           slims_only=False, aspect=None,
                           prob=stats.Binomial(), use_fdr=True,
                           progress_callback=None):
        """ Return a dictionary of enriched terms, with tuples of
        (list_of_genes, p_value, reference_count) for items and term
        ids as keys. P-Values are FDR adjusted if use_fdr is True (default).

        :param genes: List of genes
        :param reference:
            List of genes (if None all genes included in the annotations
            will be used).
        :param evidence_codes: List of evidence codes to consider.
        :param slims_only: If `True` return only slim terms.
        :param aspect:
            Which aspects to use. Use all by default. "P", "F", "C"
            or a set containing these elements.

        """
        revGenesDict = self.get_gene_names_translator(genes)
        genes = set(revGenesDict.keys())
        if reference:
            refGenesDict = self.get_gene_names_translator(reference)
            reference = set(refGenesDict.keys())
        else:
            reference = self.gene_names

        if aspect == None:
            aspects_set = set(["P", "C", "F"])
        elif isinstance(aspect, basestring):
            aspects_set = set([aspect])
        else:
            aspects_set = aspect

        evidence_codes = set(evidence_codes or evidenceDict.keys())
        annotations = [ann
                       for gene in genes for ann in self.gene_annotations[gene]
                       if ann.Evidence_Code in evidence_codes and
                       ann.Aspect in aspects_set]

        refAnnotations = set(
            [ann
             for gene in reference for ann in self.gene_annotations[gene]
             if ann.Evidence_Code in evidence_codes and
             ann.Aspect in aspects_set]
        )

        annotationsDict = defaultdict(set)
        for ann in annotations:
            annotationsDict[ann.GO_ID].add(ann)

        self._ensure_ontology()
        if slims_only and not self.ontology.slimsSubset:
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

        milestones = progress_bar_milestones(len(terms), 100)
        for i, term in enumerate(terms):
            if slims_only and term not in self.ontology.slimsSubset:
                continue
            allAnnotations = self.get_all_annotations(term).intersection(refAnnotations)
##            allAnnotations.intersection_update(refAnnotations)
            allAnnotatedGenes = set([ann.geneName for ann in allAnnotations])
            mappedGenes = genes.intersection(allAnnotatedGenes)

            if len(reference) > len(allAnnotatedGenes):
                mappedReferenceGenes = reference.intersection(allAnnotatedGenes)
            else:
                mappedReferenceGenes = allAnnotatedGenes.intersection(reference)
            res[term] = ([revGenesDict[g] for g in mappedGenes],
                         prob.p_value(len(mappedGenes), len(reference),
                                      len(mappedReferenceGenes), len(genes)),
                         len(mappedReferenceGenes))
            if progress_callback and i in milestones:
                progress_callback(100.0 * i / len(terms))
        if use_fdr:
            res = sorted(res.items(), key=lambda (_1, (_2, p, _3)): p)
            res = dict([(id, (genes, p, ref))
                        for (id, (genes, _, ref)), p in
                        zip(res, stats.FDR([p for _, (_, p, _) in res]))])
        return res

    @deprecated_keywords(
        {"directAnnotationOnly": "direct_annotation_only",
         "evidenceCodes": "evidence_codes",
         "progressCallback": "progress_callback"})
    def get_annotated_terms(self, genes, direct_annotation_only=False,
                            evidence_codes=None, progress_callback=None):
        """Return all terms that are annotated by genes with evidence_codes.
        """
        genes = [genes] if type(genes) == str else genes
        revGenesDict = self.get_gene_names_translator(genes)
        genes = set(revGenesDict.keys())
        evidence_codes = set(evidence_codes or evidenceDict.keys())
        annotations = [ann for gene in genes for ann in self.gene_annotations[gene]
                       if ann.Evidence_Code in evidence_codes]
        dd = defaultdict(set)
        for ann in annotations:
            dd[ann.GO_ID].add(revGenesDict.get(ann.geneName, ann.geneName))
        if not direct_annotation_only:
            self._ensure_ontology()
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
                termAnnots = self.get_all_annotations(term).intersection(annotations)
##                termAnnots.intersection_update(annotations)
                dd[term].update([revGenesDict.get(ann.geneName, ann.geneName)
                                 for ann in termAnnots])
        return dict(dd)

    @deprecated_keywords(
        {"clusterSize": "cluster_size", "refSize": "ref_size"})
    def draw_enrichment_graph(self, terms, cluster_size, ref_size=None,
                              file="graph.png", width=None, height=None,
                              precison=3):
        ref_size = len(self.gene_names) if ref_size == None else ref_size
        sortedterms = sorted(terms.items(), key=lambda term: term[1][1])
        fdr = dict(zip([t[0] for t in sortedterms],
                       stats.FDR([t[1][1] for t in sortedterms])))
        termsList = [(term,
                      ((float(len(terms[term][0])) / cluster_size) /
                       (float(terms[term][2]) / ref_size)),
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
        self.add_annotation(line)

    def append(self, line):
        """ Add one annotation
        """
        self.add_annotation(line)

    def extend(self, lines):
        """ Add multiple annotations
        """
        for line in lines:
            self.add_annotation(line)

    def remap_genes(self, map):
        for gene in map:
            annotations = self.gene_annotations[gene]
            for ann in annotations:
                for name in map[gene]:
                    ann1 = ann._replace(DB_Object_Symbol=name)
                    self.add(ann1)
        self.genematcher = obiGene.GMDirect()
        self._gene_names = None
        self.genematcher.set_targets(self.gene_names)

    RemapGenes = remap_genes

    @staticmethod
    @deprecated_keywords({"progressCallback": "progress_callback"})
    def download_annotations(org, file, progress_callback=None):
        if isinstance(file, basestring):
            tFile = tarfile.open(file, "w:gz")
        else:
            tFile = file

        tmpDir = os.path.join(environ.buffer_dir, "tmp_go/")
        try:
            os.mkdir(tmpDir)
        except Exception:
            pass
        fileName = "gene_association." + org + ".gz"

        from Orange.utils import wget
        wget("http://www.geneontology.org/gene-associations/" + fileName,
             directory=tmpDir,
             progress=progress_callback)

        gzFile = GzipFile(os.path.join(tmpDir, fileName), "r")
        file = open(os.path.join(tmpDir, "gene_association." + org), "w")
        file.writelines(gzFile.readlines())
        file.flush()
        file.close()

        tFile.add(os.path.join(tmpDir, "gene_association." + org),
                  "gene_association")
        annotation = Annotations(os.path.join(tmpDir, "gene_association." + org),
                    genematcher=obiGene.GMDirect(), progress_callback=progress_callback)
        cPickle.dump(annotation.gene_names, open(os.path.join(tmpDir, "gene_names.pickle"), "wb"))
        tFile.add(os.path.join(tmpDir, "gene_names.pickle"), "gene_names.pickle")
        tFile.close()
        os.remove(os.path.join(tmpDir, "gene_association." + org))
        os.remove(os.path.join(tmpDir, "gene_names.pickle"))

    DownloadAnnotations = download_annotations

    @staticmethod
    @deprecated_keywords({"progressCallback": "progress_callback"})
    def download_annotations_at_rev(org, rev, filename=None,
                                    progress_callback=None):
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

    DownloadAnnotationsAtRev = download_annotations_at_rev

from .taxonomy import pickled_cache


@pickled_cache(None, [("GO", "taxonomy.pickle"),
                      ("Taxonomy", "ncbi_taxonomy.tar.gz")])
def organism_name_search(name):
    return Annotations.organism_name_search(name)


@deprecated_keywords({"maxPValue": "p_value"})
def filter_by_p_value(terms, p_value=0.01):
    """ Filters the terms by the p-value. Assumes terms is a dict with
    the same structure as returned from get_enriched_terms.

    """
    return dict(filter(lambda (k, e): e[1] <= p_value, terms.items()))

filterByPValue = filter_by_p_value


@deprecated_keywords({"minF": "min_freq"})
def filter_by_frequency(terms, min_freq=2):
    """ Filters the terms by the cluster frequency. Asumes terms is
    a dict with the same structure as returned from get_enriched_terms.

    """
    return dict(filter(lambda (k, e): len(e[0]) >= min_freq, terms.items()))

filterByFrequency = filter_by_frequency


@deprecated_keywords({"minF": "min_freq"})
def filter_by_ref_frequency(terms, min_freq=4):
    """ Filters the terms by the reference frequency. Assumes terms is
    a dict with the same structure as returned from get_enriched_terms.

    """
    return dict(filter(lambda (k, e): e[2] >= min_freq, terms.items()))

filterByRefFrequency = filter_by_ref_frequency


def draw_enrichment_graph(enriched, file="graph.png", width=None, height=None,
                          header=None, ontology=None, precison=3):
    file = open(file, "wb") if type(file) == str else file
    _draw_enrichment_graph_tostream(enriched, file, width, height, header,
                                    ontology, precison)

drawEnrichmentGraph = draw_enrichment_graph


def _draw_enrichment_graph_tostream(enriched, fh, width, height, header=None,
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

    _draw_enrichment_graph_PIL_tostream(termsList, header, fh, width, height)

drawEnrichmentGraph_tostreamMk2 = _draw_enrichment_graph_tostream


def _draw_enrichment_graph_PIL_tostream(termsList, headers, fh, width=None,
                                        height=None):
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

drawEnrichmentGraphPIL_tostream = _draw_enrichment_graph_PIL_tostream


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
    """Maps NCBI taxonomy ids to corresponding GO organism codes.
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
            path = serverfiles.localpath_download("GO", "taxonomy.pickle")
            if os.path.isfile(path):
                self.tax = cPickle.load(open(path, "rb"))
            else:
                serverfiles.download("GO", "taxonomy.pickle")
                self.tax = cPickle.load(open(path, "rb"))

    def __getitem__(self, key):
        key = self.common_org_map.get(key, key)
        return self.code_map[key]

    def keys(self):
        return list(set(self.common_org_map.keys() + self.code_map.keys()))


def from_taxid(id):
    """ Return a set of GO organism codes that correspond to NCBI taxonomy id.
    """
    return Taxonomy()[id]


def to_taxid(db_code):
    """ Return a set of NCBI taxonomy ids from db_code GO organism annotations.
    """
    r = [key for key, val in Taxonomy().code_map.items() if db_code == val]
    return set(r)


def _test2():
    o = Ontology()
    a = Annotations("human", ontology=o)
    clusterGenes = sorted(a.gene_names)[:100]
    for i in range(10):
        genes = clusterGenes[i * 10: (i + 1) * 10]
        a.get_enriched_terms(genes, aspect=["P"])
        a.get_enriched_terms(genes, aspect=["C"])
        a.get_enriched_terms(genes, aspect=["F"])
        print i

    terms = a.get_enriched_terms(clusterGenes, aspect=["P"])
    a.get_annotated_terms(clusterGenes)

    a.draw_enrichment_graph(filterByPValue(terms), len(clusterGenes), len(a.gene_names))


def _test3():
    o = Ontology()
    a = Annotations("sgd", ontology=o)
    clusterGenes = sorted(a.gene_names)[:1] + sorted(a.gene_names)[-1:]
    exonMap = dict([(gene, [gene + "_E%i" % i for i in range(10)]) for gene in a.gene_names])
    a.RemapGenes(exonMap)

    terms = a.get_enriched_terms(exonMap.values()[0][:2] + exonMap.values()[-1][2:])

    print terms

    a.draw_enrichment_graph(filterByPValue(terms, maxPValue=0.1), len(clusterGenes), len(a.gene_names))

if __name__ == "__main__":
    _test2()
