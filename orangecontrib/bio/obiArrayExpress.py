from __future__ import absolute_import

from .arrayexpress import *

import warnings

"""
Gene Expression Atlas
---------------------

.. WARNING:: Deprecated, use ``obiGeneAtlas`` instead.

`Gene Expression Atlas <http://www.ebi.ac.uk/gxa/>`_ is a curated subset of
gene expression experiments in Array Express Archive.

Use `query_atlas_simple` for simple querys.

Example (query human genes for experiments in which they are up regulated) ::

    >>> obiArrayExpress.query_atlas_simple(genes=["SORL1", "PSIP1", "CDKN1C"], regulation="up", organism="Homo sapiens")
    {u'...

Or use the `AtlasCondition` subclasses in this module to construct a more
advanced query and use the `query_atlas` function.

Example (query human genes annotated to the GO term 'transporter activity'
that are up regulated in the liver in at least three experiments) ::

    >>> go_cond = AtlasConditionGeneProperty("Goterm", "Is", "transporter activity")
    >>> liver_cond = AtlasConditionExperimentalFactor("Organism_part", "up", 3, "liver")
    >>> org_cond = AtlasConditionOrganism("Homo sapiens")
    >>> cond_list = AtlasConditionList([go_cond, liver_cond, org_cond])
    >>> query_atlas(cond_list)
    {u'...

"""


def _open_shelve(filename, flag="r"):
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    exists = os.path.exists(filename)
    if flag in ["r", "w"] and not exists:
        # needs to be created first
        # XXX: Race condition
        s = shelve.open(filename, "c")
        s.close()

    return shelve.open(filename, flag)


@contextmanager
def _fake_closing(obj):
    yield obj


class GeneExpressionAtlasConenction(object):

    """
    A connection to Gene Expression Atlas database.

    :param address:
        Address of the GXA server (default: http://www.ebi.ac.uk/gxa/api/deprecated).
    :param timeout:
        Socket timeout (default 30).
    :param cache:
        A dict like object to use as a cache.

    """
    DEFAULT_ADDRESS = "http://www-test.ebi.ac.uk/gxa/api/deprecated"

    DEFAULT_CACHE = serverfiles.localpath(
        "ArrayExpress", "GeneAtlasCache.shelve")

    def __init__(self, address=None, timeout=30, cache=None):
        """
        Initialize the connection.

        """
        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.timeout = timeout
        if cache is None:
            cache = self.DEFAULT_CACHE

        self.cache = cache

    def query(self, condition, format="json", start=None, rows=None, indent=False):
        url = self.address + "?" + condition.rest()
        if start is not None and rows is not None:
            url += "&start={0}&rows={1}".format(start, rows)
        url += "&format={0}".format(format)
        if indent:
            url += "&indent"
#        print url
        if self.cache is not None:
            return self._query_cached(url)
        else:
            return urllib2.urlopen(url)

    def _query_cached(self, url):
        if self.cache is not None:
            with self.open_cache("r") as cache:
                if url in cache:
                    return StringIO(cache[url])

            response = urllib2.urlopen(url)
            contents = response.read()
            with self.open_cache("w") as cache:
                cache[url] = contents

            return StringIO(contents)
        else:
            return urllib2.urlopen(url)

    def open_cache(self, flag="r"):
        """
        Return a context manager for a dict like object.
        """
        if isinstance(self.cache, basestring):
            try:
                return closing(_open_shelve(self.cache, flag))
            except Exception:
                return _fake_closing({})
        else:
            return _fake_closing(self.cache)


# Names of all Gene Property filter names
GENE_FILTERS = \
    ["Name",  # Gene name
     "Goterm",  # Gene Ontology Term
     "Interproterm",  # InterPro Term
     "Disease",  # Gene-Disease Assocation
     "Keyword",  # Gene Keyword
     "Protein",  # Protein

     "Dbxref",  # Other Database Cross-Refs
     "Embl",  # EMBL-Bank ID
     "Ensfamily",  # Ensembl Family
     "Ensgene",  # Ensembl Gene ID

     "Ensprotein",  # Ensembl Protein ID
     "Enstranscript",  # Ensembl Transcript ID
     "Goid",  # Gene Ontology ID
     "Image",  # IMAGE ID
     "Interproid",  # InterPro ID
     "Locuslink",  # Entrez Gene ID

     "Omimid",  # OMIM ID
     "Orf",  # ORF
     "Refseq",  # RefSeq ID
     "Unigene",  # UniGene ID
     "Uniprot",  # UniProt Accession

     "Hmdb",  # HMDB ID
     "Chebi",  # ChEBI ID
     "Cas",  # CAS
     "Uniprotmetenz",  # Uniprotmetenz
     "Gene",  # Gene Name or Identifier
     "Synonym",  # Gene Synonym
     ]

# Valid Gene Property filter qualifiers
GENE_FILTER_QUALIFIERS =\
    ["Is",
     "IsNot"
     ]

# Organisms in the Atlas
ATLAS_ORGANISMS = \
    ["Anopheles gambiae",
     "Arabidopsis thaliana",
     "Bos taurus",
     "Caenorhabditis elegans",
     "Danio rerio",
     "Drosophila melanogaster",
     "Epstein barr virus",
     "Gallus gallus",
     "Homo sapiens",
     "Human cytomegalovirus",
     "Kaposi sarcoma-associated herpesvirus",
     "Mus musculus",
     "Rattus norvegicus",
     "Saccharomyces cerevisiae",
     "Schizosaccharomyces pombe",
     #     "Unknown",
     "Xenopus laevis"
     ]

#_COMMON_TAXIDS = \
#    {"Anopheles gambiae",
#     "Arabidopsis thaliana",
#     "Bos taurus",
#     "Caenorhabditis elegans",
#     "Danio rerio",
#     "Drosophila melanogaster",
#     "Epstein barr virus",
#     "Gallus gallus",
#     "Homo sapiens",
#     "Human cytomegalovirus",
#     "Kaposi sarcoma-associated herpesvirus",
#     "Mus musculus",
#     "Rattus norvegicus",
#     "Saccharomyces cerevisiae",
#     "Schizosaccharomyces pombe",
# "Unknown",
#     "Xenopus laevis"
#     }


def ef_ontology():
    """ Return the `EF <http://www.ebi.ac.uk/efo/>`_ (Experimental Factor) ontology
    """
    from . import obiOntology
#    return obiOntology.OBOOntology(urllib2.urlopen("http://efo.svn.sourceforge.net/svnroot/efo/trunk/src/efoinobo/efo.obo"))
    # Should this be in the OBOFoundry (Ontology) domain
    try:
        file = open(serverfiles.localpath_download("ArrayExpress", "efo.obo"), "rb")
    except urllib2.HTTPError:
        file = urllib2.urlopen("http://efo.svn.sourceforge.net/svnroot/efo/trunk/src/efoinobo/efo.obo")
    return obiOntology.OBOOntology(file)


class AtlasCondition(object):

    """ Base class for Gene Expression Atlas query condition
    """

    def validate(self):
        """ Validate condition in a subclass.
        """
        raise NotImplementedError

    def rest(self):
        """ Return a REST query part in a subclass.
        """
        raise NotImplementedError


class AtlasConditionList(list, AtlasCondition):

    """ A list of AtlasCondition instances.
    """

    def validate(self):
        for item in self:
            item.validate()

    def rest(self):
        return "&".join(cond.rest() for cond in self)


class AtlasConditionGeneProperty(AtlasCondition):

    """ An atlas gene filter condition.

    :param property: Property of the gene. If None or "" all properties 
        will be searched.
    :param qualifier: Qualifier can be 'Is' or 'IsNot'
    :param value: The value to search for.

    Example ::

        >>> # Condition on a gene name
        >>> condition = AtlasConditionGeneProperty("Name", "Is", "AS3MT")
        >>> # Condition on genes from a GO Term
        >>> condition = AtlasConditionGeneProperty("Goterm", "Is", "p53 binding")
        >>> # Condition on disease association
        >>> condition = AtlasConditionGeneProperty("Disease", "Is", "cancer")

    """

    def __init__(self, property, qualifier, value):
        self.property = property or ""
        self.qualifier = qualifier
        if isinstance(value, basestring):
            self.value = value.replace(" ", "+")
        elif isinstance(value, list):
            self.value = "+".join(value)
        else:
            raise ValueError(value)

        self.validate()

    def validate(self):
        assert(self.property in GENE_FILTERS + [""])
        assert(self.qualifier in GENE_FILTER_QUALIFIERS + [""])

    def rest(self):
        return "gene{property}{qualifier}={value}".format(**self.__dict__)


class AtlasConditionExperimentalFactor(AtlasCondition):

    """ An atlas experimental factor filter condition.

    :param factor: EFO experiamntal factor
    :param regulation: "up", "down", "updown", "any" or "none"
    :param n: Minimum number of of experimants with this condition
    :param value: Experimantal factor value

    Example ::

        >>> # Any genes up regulated in at least 3 experiments involving cancer.
        >>> condition = AtlasConditionExperimentalFactor("", "up", 3, "cancer")
        >>> # Only genes which are up/down regulated in the heart in at least one experiment. 
        >>> condition = AtlasConditionExperimentalFactor("Organism_part", "updown", 1, "heart")

    """

    def __init__(self, factor, regulation, n, value):
        self.factor = factor
        self.regulation = regulation
        self.n = n
        self.value = value
        self.validate()

    def validate(self):
        # TODO: validate the factor and value
#        assert(self.factor in ef_ontology())
        assert(self.regulation in ["up", "down", "updown"])

    def rest(self):
        return "{regulation}{n}In{factor}={value}".format(**self.__dict__)


class AtlasConditionOrganism(AtlasCondition):

    """ Condition on organism.
    """

    def __init__(self, organism):
        self.organism = organism
        self.validate()

    def validate(self):
        assert(self.organism in ATLAS_ORGANISMS)

    def rest(self):
        return "species={0}".format(self.organism.replace(" ", "+").lower())


class AtlasConditionExperiment(AtlasCondition):

    """ Condition on experiement

    :param property: Property of the experiment. If None or "" all properties
        will be searched.
    :param qualifier: Qualifier can be 'Has' or 'HasNot'
    :param value: The value to search for.

    Example ::

        >>> # Condition on a experiemnt acession
        >>> condition = AtlasConditionExperiment("", "", "E-GEOD-24283")
        >>> # Condition on experiments involving lung
        >>> condition = AtlasConditionExperiment("Organism_part", "Has", "lung")

    """
#    EXPERIMENT_FILTERS = [
#                "Organism"
#                "Factor"]

    EXPERIMENT_FILTER_QUALIFIERS = [
        "Has",
        "HasNot"]

    def __init__(self, property, qualifier, value):
        self.property = property
        self.qualifier = qualifier
        if isinstance(value, basestring):
            self.value = value.replace(" ", "+")
        elif isinstance(value, list):
            self.value = "+".join(value)
        else:
            raise ValueError(value)

        self.validate()

    def validate(self):
        # TODO: check to EFO factors
#        assert(self.property in EXPERIMENT_FILTERS + [""])
        assert(self.qualifier in self.EXPERIMENT_FILTER_QUALIFIERS + [""])

    def rest(self):
        return "experiment{property}{qualifier}={value}".format(**self.__dict__)


class GeneAtlasError(ValueError):

    """ An error response from the Atlas server.
    """
    pass


def __check_atlas_error_json(response):
    if "error" in response:
        raise GeneAtlasError(response["error"])
    return response


def __check_atlas_error_xml(response):
    error = response.find("error")
    if error is not None:
        raise GeneAtlasError(error.text)
    return response


def query_atlas_simple(genes=None, regulation=None, organism=None,
                       condition=None, format="json", start=None,
                       rows=None):
    """ A simple Atlas query.

    :param genes: A list of gene names to search for.
    :param regulation: Search for experiments in which `genes` are "up",
        "down", "updown" or "none" regulated. If None all experiments
        are searched.
    :param organism: Search experiments for organism. If None all experiments
        are searched.
    :param condition: An EFO factor value (e.g. "brain")

    Example ::

        >>> query_atlas_simple(genes=['Pou5f1', 'Dppa3'], organism="Mus musculus")
        {u'...

        >>> query_atlas_simple(genes=['Pou5f1', 'Dppa3'], regulation="up", organism="Mus musculus")
        {u'...

        >>> query_atlas_simple(genes=['Pou5f1', 'Dppa3'], regulation="up", condition="liver", organism="Mus musculus")
        {u'...

    """
    warnings.warn("Use 'obiGeneAtlas.run_simple_query' instead.", DeprecationWarning)
    conditions = AtlasConditionList()
    if genes:
        conditions.append(AtlasConditionGeneProperty("Gene", "Is", genes))
    if regulation or condition:
        regulation = "any" if regulation is None else regulation
        condition = "" if condition is None else condition
        conditions.append(AtlasConditionExperimentalFactor("", regulation, 1, condition))
    if organism:
        conditions.append(AtlasConditionOrganism(organism))

    connection = GeneExpressionAtlasConenction()
    results = connection.query(conditions, format=format, start=start,
                               rows=rows)
    if format == "json":
        return parse_json(results)
    else:
        return parse_xml(results)


# TODO: can this be implemented query_atlas(organism="...", Locuslink="...", Chebi="...", up3InCompound="..." downInEFO="...")
# Need a full list of accepted factors


def query_atlas(condition, format="json", start=None, rows=None, indent=False, connection=None):
    """ Query Atlas based on a `condition` (instance of AtlasCondition)

    Example ::

        >>> condition1 = AtlasConditionGeneProperty("Goterm", "Is", "p53 binding")
        >>> condition2 = AtlasConditionExperimentalFactor("Organism_part", "up", 3, "heart")
        >>> condition = AtlasConditionList([condition1, condition2])
        >>> query_atlas(condition)
        {u'...

    """
    warnings.warn("Use 'obiGeneAtlas.run_query' instead.", DeprecationWarning)
    if connection is None:
        connection = GeneExpressionAtlasConenction()
    results = connection.query(condition, format=format, start=start,
                               rows=rows, indent=indent)
    if format == "json":
        response = parse_json(results)
        return __check_atlas_error_json(response)
    else:
        response = parse_xml(results)
        return __check_atlas_error_xml(response)


def get_atlas_summary(genes, organism, connection=None):
    """ Return 3 dictionaries containing a summary of atlas information
    about three experimental factors:

        - Organism Part (OP) 
        - Disease State (DS)
        - Cell type (CT)

    Each dictionary contains query genes as keys. Values are dictionaries
    mapping factor values to a 2-tuple containig the count of up regulated
    and down regulated experiments.

    Example ::

        >>> get_atlas_summary(["RUNX1"], "Homo sapiens")
        ({u'RUNX1': ...

    """
    warnings.warn("Use 'obiGeneAtlas.get_atlas_summary' instead.", DeprecationWarning)
    genes_condition = AtlasConditionGeneProperty("Gene", "Is", genes)
    org_condition = AtlasConditionOrganism(organism)
    condition = AtlasConditionList([genes_condition, org_condition])
    result = query_atlas(condition, format="json", connection=connection)

    org_part = collect_ef_summary(result, "organism_part")
    disease_state = collect_ef_summary(result, "disease_state")
    cell_type = collect_ef_summary(result, "cell_type")

    return dict(org_part), dict(disease_state), dict(cell_type)


def collect_ef_summary(info, ef, summary=None):
    """ Collect the results summary from query_atlas, result for experimental
    factor `ef`.
    """
    if summary is None:
        summary = defaultdict(dict)

    results = info["results"]
    for res in results:
        gene = res["gene"]
        expressions = res["expressions"]
        for expression in expressions:
            if expression["ef"] == ef:
                efv = expression["efv"]
                updown = (expression["upExperiments"],
                          expression["downExperiments"]
                          )

                if any(updown):
                    summary[gene["name"]][efv] = updown

    return summary
