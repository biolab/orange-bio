"""
========================================
Gene Expression Atlas (``obiGeneAtlas``)
========================================

.. warning::

    Deprecated!!!
    Gene Expression Atlas REST api has been deprecated and will be removed
    in the future.


Interface to Gene Expression Atlas.

`Gene Expression Atlas <http://www.ebi.ac.uk/gxa/>`_ is a curated subset of
gene expression experiments in Array Express Archive.


.. autofunction:: gene_expression_atlas

.. autofunction:: default_gene_matcher

.. autofunction:: to_taxid

"""

from __future__ import absolute_import

import os
import shelve
import warnings
from collections import defaultdict, namedtuple
from contextlib import closing, contextmanager

from Orange.utils import serverfiles

from . import obiGene

GeneResults = namedtuple("GeneResults", "id name synonyms expressions")
ExpressionResults = namedtuple("ExpressionResults", "ef efv up down experiments")
ExperimentExpression = namedtuple("ExperimentExpression", "accession expression pvalue")

##
GeneAtlasResult = GeneResults
AtlasExpressions = ExpressionResults
AtlasExperiment = ExperimentExpression
##

CACHE_VERSION = 1


def _cache(name="AtlasGeneResult.shelve"):
    """ Return a open cache instance (a shelve object).
    """
    if not os.path.exists(serverfiles.localpath("GeneAtlas")):
        try:
            os.makedirs(serverfiles.localpath("GeneAtlas"))
        except OSError:
            pass
    cache = shelve.open(serverfiles.localpath("GeneAtlas", name))
    if cache.get(name + "__CACHE_VERSION__", None) == CACHE_VERSION:
        return cache
    else:
        cache.close()
        cache = shelve.open(serverfiles.localpath("GeneAtlas", name), "n")
        cache[name + "__CACHE_VERSION__"] = CACHE_VERSION
        return cache


SLEEP_TIME_MULTIPLIER = 3.0

def gene_expression_atlas(genes, progress_callback=None):
    """ Return GeneResults instances for genes (genes must be valid ensembl ids).
    """
    import time
    genes = list(genes)
    result_dict = {}
    genes_not_cached = []
    # See which genes are already cached
    with closing(_cache()) as cache:
        for gene in genes:
            if str(gene) in cache:
                result_dict[gene] = cache[str(gene)]
            else:
                genes_not_cached.append(gene)
    
    batch_size = 10
    start = 0
    res = []
    while start < len(genes_not_cached):
        batch = genes_not_cached[start: start + batch_size]
        start += batch_size
        start_time = time.time()
        batch_res = batch_gene_atlas_expression(batch)
        # Cache the new results.
        # TODO: handle genes without any results.
        genes_with_no_results = set(batch) - set(r.id for r in batch_res) 
        with closing(_cache()) as cache:
            for atlas_res in batch_res:
                cache[str(atlas_res.id)] = atlas_res
                result_dict[atlas_res.id] = atlas_res
            for g in genes_with_no_results:
                cache[str(g)] = None
        res.extend(batch_res)
        # Sleep
        if start % (batch_size * 10) == 0:
            # every 10 batches wait one minute before continuing. 
            time.sleep(60)
        else:
            time.sleep(min(20.0, SLEEP_TIME_MULTIPLIER*(time.time() - start_time)))
            
        if progress_callback:
            progress_callback(100.0 * start / len(genes_not_cached))
    
    return [result_dict.get(g, None) for g in genes]

    
def batch_gene_atlas_expression(genes):
    cond = GenePropertyCondition("Ensgene", "Is", genes)
    res = run_query(cond, format="json")
    results = res["results"]
    results_genes = []
    for one_result in results:
        gene = one_result["gene"]
        id = gene["id"]
        name = gene["name"]
        synonyms = gene.get("synonyms", [])
        expressions = one_result["expressions"]
        result_expressions = []
        for expression in expressions:
            ef = expression["ef"]
            efv = expression["efv"]
            up = expression["upExperiments"]
            down = expression["downExperiments"]
            experiments = expression["experiments"]
            result_experiment = []
            for exp in experiments:
                if "accession" in exp:
                    exp_accession = exp["accession"]
                elif "experimentAccession" in exp:
                    exp_accession = exp["experimentAccession"]
                else:
                    raise KeyError()
                if "expression" in exp:
                    updown = exp["expression"]
                elif "updn" in exp:
                    updown = exp["updn"]
                else:
                    raise KeyError
                pval = exp["pvalue"]
                result_experiment.append(ExperimentExpression(exp_accession, updown, pval))
            result_expressions.append(ExpressionResults(ef, efv, up, down, result_experiment))
        results_genes.append(GeneResults(id, name, synonyms, result_expressions))
    return results_genes


def default_gene_matcher(organism):
    """ Return a default gene matcher for organism
    (targeting Ensembl gene ids).
    
    """
    taxid = to_taxid(organism)
    matcher = obiGene.matcher(
      [obiGene.GMEnsembl(taxid),
       [obiGene.GMEnsembl(taxid),
        obiGene.GMNCBI(taxid)]]
    )
    matcher.set_targets(obiGene.EnsembleGeneInfo(taxid).keys())
    return matcher

from Orange.utils import lru_cache


@lru_cache(maxsize=3)
def _cached_default_gene_matcher(organism):
    return default_gene_matcher(organism)
    

def get_atlas_summary(genes, organism, gene_matcher=None,
                      progress_callback=None):
    """ Return 3 dictionaries containing a summary of atlas information
    about three experimental factors:
    
        - Organism Part (OP) 
        - Disease State (DS)
        - Cell type (CT)
    
    Each dictionary contains query genes as keys. Values are dictionaries
    mapping factor values to a 2-tuple containig the count of up regulated
    and down regulated experiments.
    
    Example ::
    
        >>> get_atlas_summary(["ENSG00000159216"], "Homo sapiens")
        ({u'RUNX1': ...
        
    """
    if gene_matcher is None:
        gene_matcher = _cached_default_gene_matcher(organism)
        
    matched, unmatched = [], []
    for gene, match in zip(genes, map(gene_matcher.umatch, genes)):
        if match:
            matched.append(match)
        else:
            unmatched.append(gene)
    if unmatched:
        warnings.warn("Unmatched genes " + "," .join(["%r" % g for g in unmatched]))
    
    results = gene_expression_atlas(matched, progress_callback=progress_callback)
    
    def collect_ef_summary(result, ef, summary):
        for exp in result.expressions:
            if exp.ef == ef:
                if any([exp.up, exp.down]):
                    summary[result.name][exp.efv] = (exp.up, exp.down)
        
            
    op, ds, ct = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    for res in results:
        if res:
            collect_ef_summary(res, "organism_part", op)
            collect_ef_summary(res, "disease_state", ds)
            collect_ef_summary(res, "cell_type", ct)
        
    return dict(op), dict(ds), dict(ct)

def drop_none(iter):
    """ Drop all ``None`` from the iterator. 
    """
    for e in iter:
        if e is not None:
            yield e
            
def construct_atlas_gene_sets(genes, organism, factors=["organism_part",
                                    "disease_state", "cell_type"],
                              max_pvalue=1e-5):
    """ Construct gene sets for atlas experimental factor values in
    ``factors``.
    """
    results = gene_expression_atlas(genes)
    sets = defaultdict(list)
    
    for res in drop_none(results):
        for exp in res.expressions:
            if exp.ef not in factors:
                continue
            diff_exp = [e for e in exp.experiments \
                        if e.pvalue <= max_pvalue]
            if diff_exp:
                sets[exp.ef, exp.efv].append(res.id)

    organism = "+".join(organism.lower().split())
    from .obiGeneSets import GeneSets, GeneSet
    
    def display_string(name):
        return name.capitalize().replace("_", " ")
    
    gene_sets = []
    for (ef, efv), genes in sets.items():
        ef_display = display_string(ef)
        gs = GeneSet(genes, "Diff. expressed in %s=%s." % (ef_display, efv), id=ef + ":" + efv,
                     description="Diff. expressed in %s=%s" % (ef_display, efv),
                     link="http://www.ebi.ac.uk/gxa/qrs?specie_0={organism}&gprop_0=&gnot_0=&gval_0=&fact_1=&fexp_1=UPDOWN&fmex_1=&fval_1=%22{efv}%22+&view=hm".format( \
                            organism=organism, efv="+".join(efv.lower().split())),
                     hierarchy=("Gene expression atlas", ef_display))
        gene_sets.append(gs)
    return GeneSets(gene_sets)


# Mapping for common taxids from obiTaxonomy
TAXID_TO_ORG = {"": "Anopheles gambiae",
                "3702": "Arabidopsis thaliana",
                "9913": "Bos taurus",
                "6239": "Caenorhabditis elegans",
                "7955": "Danio rerio",
                "7227": "Drosophila melanogaster",
                "": "Epstein barr virus",
                "": "Gallus gallus",
                "9606": "Homo sapiens",
                "": "Human cytomegalovirus",
                "": "Kaposi sarcoma-associated herpesvirus",
                "10090": "Mus musculus",
                "10116": "Rattus norvegicus",
                "4932": "Saccharomyces cerevisiae",
                "4896": "Schizosaccharomyces pombe",
                "8355": "Xenopus laevis"
     }

def to_taxid(name):
    dd = dict((v, k) for k, v in TAXID_TO_ORG.items())
    if name in dd:
        return dd[name]
    else:
        from . import obiTaxonomy as tax
        ids = tax.to_taxid(name, mapTo=TAXID_TO_ORG.keys())
        if ids:
            return ids.pop()
        else:
            raise ValueError("Unknown organism.")


__doc__ += """\
Low level REST query interface
------------------------------

Use `query_atlas_simple` for simple querys.

Example (query human genes for experiments in which they are up regulated) ::

    >>> run_simple_query(genes=["SORL1", "PSIP1", "CDKN1C"], regulation="up", organism="Homo sapiens")
    {u'...
    
Or use the `AtlasCondition` subclasses in this module to construct a more
advanced query and use the `run_query` function.

Example (query human genes annotated to the GO term 'transporter activity'
that are up regulated in the liver in at least three experiments) ::

    >>> go_cond = GenePropertyCondition("Goterm", "Is", "transporter activity")
    >>> liver_cond = ExperimentalFactorCondition("Organism_part", "up", 3, "liver")
    >>> org_cond = OrganismCondition("Homo sapiens")
    >>> cond_list = ConditionList([go_cond, liver_cond, org_cond])
    >>> run_query(cond_list)
    {u'...

"""

import urllib2
from  StringIO import StringIO
import json
from xml.etree.ElementTree import ElementTree

parse_json = json.load


def parse_xml(stream):
    """ Parse an xml stream into an instance of xml.etree.ElementTree.ElementTree.
    """
    return ElementTree(file=stream)


class GeneExpressionAtlasConenction(object):
    """
    A connection to Gene Expression Atlas database.

    :param address:
        Address of the GXA server (default: http://www-test.ebi.ac.uk/gxa/api/deprecated).
    :param timeout:
        Socket timeout (default 30).
    :param cache:
        A dict like object to use as a cache.

    """
    DEFAULT_ADDRESS = "http://www-test.ebi.ac.uk/gxa/api/deprecated"
    DEFAULT_CACHE = serverfiles.localpath(
        "GeneAtlas", "GeneAtlasConnectionCache.shelve")

    def __init__(self, address=None, timeout=30, cache=None):

        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.timeout = timeout
        self.cache = cache if cache is not None else self.DEFAULT_CACHE

    def query(self, condition, format="json", start=None, rows=None, indent=False):
        warnings.warn(
            "The Gene Expression Atlas REST api has been deprecated and " +
            "will be removed in the future.",
            UserWarning)

        url = self.address + "?" + condition.rest()
        if start is not None and rows is not None:
            url += "&start={0}&rows={1}".format(start, rows)
        url += "&format={0}".format(format)
        if indent:
            url += "&indent"
        #print url

        if self.cache is not None:
            return self._query_cached(url, format)
        else:
            return urllib2.urlopen(url)

    def _query_cached(self, url, format):
        if self.cache is not None:
            with self.open_cache("r") as cache:
                if url in cache:
                    return StringIO(cache[url])

            response = urllib2.urlopen(url)
            contents = response.read()
            # Test if the contents is a valid json or xml string (some
            # times the stream just stops in the middle, so this makes
            # sure we don't cache an invalid response
            # TODO: what about errors (e.g. 'cannot handle the
            # query in a timely fashion'
            if format == "json":
                parse_json(StringIO(contents))
            else:
                parse_xml(StringIO(contents))

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
                return fake_closing({})
        else:
            return fake_closing(self.cache)


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
def fake_closing(obj):
    yield obj
    
    
# Names of all Gene Property filter names
GENE_FILTERS = \
    ["Name", # Gene name
     "Goterm", #Gene Ontology Term
     "Interproterm", #InterPro Term
     "Disease", #Gene-Disease Assocation
     "Keyword", #Gene Keyword
     "Protein", #Protein

     "Dbxref", #Other Database Cross-Refs
     "Embl", #EMBL-Bank ID
     "Ensfamily", #Ensembl Family
     "Ensgene", #Ensembl Gene ID

     "Ensprotein", #Ensembl Protein ID
     "Enstranscript", #Ensembl Transcript ID
     "Goid", #Gene Ontology ID
     "Image", #IMAGE ID
     "Interproid", #InterPro ID
     "Locuslink", #Entrez Gene ID

     "Omimid", #OMIM ID
     "Orf", #ORF
     "Refseq", #RefSeq ID
     "Unigene", #UniGene ID
     "Uniprot", #UniProt Accession

     "Hmdb", #HMDB ID
     "Chebi", #ChEBI ID
     "Cas", #CAS
     "Uniprotmetenz", #Uniprotmetenz
     "Gene", #Gene Name or Identifier
     "Synonym", #Gene Synonym
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
     "Bacillus subtilis",
     "Bos taurus",
     "Caenorhabditis elegans",
     "Canis familiaris",
     "Ciona intestinalis",
     "Ciona savignyi",
     "Danio rerio",
     "Dasypus novemcinctus",
     "Drosophila melanogaster",
     "Epstein barr virus",
     "Equus caballus",
     "Gallus gallus",
     "Gorilla gorilla",
     "Homo sapiens",
     "Human cytomegalovirus",
     "Human immunodeficiency virus",
     "Kaposi sarcoma-associated herpesvirus",
     "Macaca mulatta",
     "Monodelphis domestica",
     "Mus musculus",
     "Ornithorhynchus anatinus",
     "Oryza sativa",
     "Pan troglodytes",
     "Pongo abelii",
     "Populus trichocarpa",
     "Rattus norvegicus",
     "Saccharomyces cerevisiae",
     "Schizosaccharomyces pombe",
     "Sus scrofa",
     "Unknown",
     "Vitis vinifera",
     "Xenopus laevis",
     "Xenopus tropicalis"
     ]
    
def ef_ontology():
    """ Return the `EF <http://www.ebi.ac.uk/efo/>`_ (Experimental Factor) ontology
    """
    from . import obiOntology
    from Orange.utils import serverfiles
    # Should this be in the OBOFoundry (Ontology) domain
    try:
        file = open(serverfiles.localpath_download("ArrayExpress", "efo.obo"), "rb")
    except urllib2.HTTPError:
        file = urllib2.urlopen("http://efo.svn.sourceforge.net/svnroot/efo/trunk/src/efoinobo/efo.obo")
    return obiOntology.OBOOntology(file)


class Condition(object):
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
    
    
class ConditionList(list, Condition):
    """ A list of AtlasCondition instances.
    """ 
    def validate(self):
        for item in self:
            item.validate()
        
    def rest(self):
        return "&".join(cond.rest() for cond in self)


class GenePropertyCondition(Condition):
    """ An atlas gene filter condition.
    
    :param property: Property of the gene. If None or "" all properties 
        will be searched.
    :param qualifier: Qualifier can be 'Is' or 'IsNot'
    :param value: The value to search for.
    
    Example ::
    
        >>> # Condition on a gene name
        >>> condition = GenePropertyCondition("Name", "Is", "AS3MT")
        >>> # Condition on genes from a GO Term
        >>> condition = GenePropertyCondition("Goterm", "Is", "p53 binding")
        >>> # Condition on disease association
        >>> condition = GenePropertyCondition("Disease", "Is", "cancer")
        
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
        
        
class ExperimentalFactorCondition(Condition):
    """ An atlas experimental factor filter condition.
    
    :param factor: EFO experiamntal factor
    :param regulation: "up", "down", "updown", "any" or "none"
    :param n: Minimum number of of experimants with this condition
    :param value: Experimantal factor value
    
    Example ::
    
        >>> # Any genes up regulated in at least 3 experiments involving cancer.
        >>> condition = ExperimentalFactorCondition("", "up", 3, "cancer")
        >>> # Only genes which are up/down regulated in the heart in at least one experiment. 
        >>> condition = ExperimentalFactorCondition("Organism_part", "updown", 1, "heart")
        
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
        assert(self.regulation in ["up", "down", "updown", "none", "any"])
        
    def rest(self):
        return "{regulation}{n}In{factor}={value}".format(**self.__dict__)
        
        
class OrganismCondition(Condition):
    """ Condition on organism.
    """
    def __init__(self, organism):
        self.organism = organism
        self.validate()
        
    def validate(self):
        assert(self.organism in ATLAS_ORGANISMS)
        
    def rest(self):
        return "species={0}".format(self.organism.replace(" ", "+").lower())
        
        
class ExperimentCondition(Condition):
    """ Condition on experiement
    
    :param property: Property of the experiment. If None or "" all properties 
        will be searched.
    :param qualifier: Qualifier can be 'Has' or 'HasNot'
    :param value: The value to search for.
    
    Example ::
    
        >>> # Condition on a experiemnt acession
        >>> condition = ExperimentCondition("", "", "E-GEOD-24283")
        >>> # Condition on experiments involving lung
        >>> condition = ExperimentCondition("Organism_part", "Has", "lung")
        
    """
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
        
        
class GeneExpressionAtlasError(Exception):
    """ An error response from the Atlas server.
    """
    pass
    
    
def __check_atlas_error_json(response):
    if "error" in response:
        raise GeneExpressionAtlasError(response["error"])
    return response
 
     
def __check_atlas_error_xml(response):
    error = response.find("error")
    if error is not None:
        raise GeneExpressionAtlasError(error.text)
    return response
    
        
def run_simple_query(genes=None, regulation=None, organism=None,
                     condition=None, format="json", start=None,
                     rows=None):
    """ A simple Gene Atlas query.
    
    :param genes: A list of gene names to search for.
    :param regulation: Search for experiments in which `genes` are "up",
        "down", "updown", "any" or "none" regulated. If "any" all experiments
        are searched.
    :param organism: Search experiments for organism. If None all experiments
        are searched.
    :param condition: An EFO factor value (e.g. "brain")
    
    Example ::
        
        >>> run_simple_query(genes=['Pou5f1', 'Dppa3'], organism="Mus musculus")
        {u'...
        
        >>> run_simple_query(genes=['Pou5f1', 'Dppa3'], regulation="up", organism="Mus musculus")
        {u'...
        
        >>> run_simple_query(genes=['Pou5f1', 'Dppa3'], regulation="up", condition="liver", organism="Mus musculus")
        {u'...
        
    """
    conditions = ConditionList()
    if genes:
        conditions.append(GenePropertyCondition("Gene", "Is", genes))
    if regulation or condition:
        regulation = "any" if regulation is None else regulation
        condition = "" if condition is None else condition
        conditions.append(ExperimentalFactorCondition("", regulation, 1, condition))
    if organism:
        conditions.append(OrganismCondition(organism))
        
    connection = GeneExpressionAtlasConenction()
    results = connection.query(conditions, format=format, start=start,
                               rows=rows)
    if format == "json":
        return parse_json(results)
    else:
        return parse_xml(results)

"""\
.. todo:: can this be implemented query_atlas(organism="...", Locuslink="...", Chebi="...", up3InCompound="..." downInEFO="...")
      Need a full list of accepted factors 
"""

def run_query(condition, format="json", start=None, rows=None, indent=False, connection=None):
    """ Query Atlas based on a `condition` (instance of :class:`Condition`)
    
    Example ::
        
        >>> condition1 = GenePropertyCondition("Goterm", "Is", "p53 binding")
        >>> condition2 = ExperimentalFactorCondition("Organism_part", "up", 3, "heart")
        >>> condition = ConditionList([condition1, condition2])
        >>> run_query(condition)
        {u'...
        
    """
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
    
def test():
    from pprint import pprint    
    pprint(get_atlas_summary(['Pou5f1', 'Dppa3'], 'Mus musculus'))
       
    pprint(get_atlas_summary(['PDLIM5', 'FGFR2' ], 'Homo sapiens'))
    import doctest 
    doctest.testmod(optionflags=doctest.ELLIPSIS)
    
if __name__ == "__main__":
    test()
