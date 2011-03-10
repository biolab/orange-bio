"""
Array Express
-------------

A python module for accessing the ArrayExpress and GeneExpressionAtlas web services.

Example::

    >>> import obiArrayExpress
    >>> obiArrayExpress.query_experiments(accession='E-MEXP-31')
    <addinfourl at ...
    
    
    
"""

import os, sys
import urllib2

import orngEnviron
import warnings
import posixpath
import shelve
from collections import defaultdict

ARRAYEXPRESS_FIELDS = \
    ["accession",
     "array",
     "ef",
     "efv",
     "expdesign",
     "exptype",
     "gxa",
     "pmid",
     "sa",
     "species",
    ]

class ArrayExpressConnection(object):
    """ A connection to the ArrayExpress database with query caching
    """
    
    try:
        DEFAULT_CACHE = shelve.open(os.path.join(orngEnviron.bufferDir, "ArrayExpressCache.shelve"))
    except ImportError:
        warnings.warn("Could not load persistent cache!")
        DEFAULT_CACHE = {}
    
    DEFAULT_ADDRESS = "http://www.ebi.ac.uk/arrayexpress/{format}/v2/"
    DEFAULT_FORMAT = "json"
    
    # Order of arguments in the query
    _ARGS_ORDER = ["keywords", "species", "array"]
    
    def __init__(self, address=None, cache=None, timeout=30):
        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.cache = cache if cache is not None else self.DEFAULT_CACHE
        self.timeout = timeout
        
    def format_query(self, **kwargs):
        def format_species(val):
            return '"%s"' % val.lower()
        # TODO: range values (e.g. [1 TO 2]), handle AND, OR, +, check for valid keys
        formaters = {"species": format_species,
                     }
        parts = []
        arg_items = kwargs.items()
        ordered = sorted(arg_items, key=lambda arg: self._ARGS_ORDER.index(arg[0]) \
                         if arg[0] in self._ARGS_ORDER else 100)
        for key, value in kwargs.iteritems():
            fmt = formaters.get("key", lambda val: val)
            value = fmt(value)
            parts.append("{0}={1}".format(key, value)) 
        return "&".join(parts)
        
    def query_url(self, what="experiments", **kwargs):
        """ Return a formated query url for the calls arguments
        
        Example::
            >>> conn.query_url(accession="E-MEXP-31")
            'http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?accession=E-MEXP-31'
            
        """
        query = self.format_query(**kwargs)
        url = posixpath.join(self.address, what)
        url = url.format(format=kwargs.get("format", self.DEFAULT_FORMAT))
        url = url + "?" + query
#        print url
        url = url.replace(" ", "%20")
        return url
    
    def query_url_experiments(self, **kwargs):
        """ Return a formated experiments query url for the calls arguments
        """
        return self.query_url("experiments", **kwargs)
    
    def query_url_files(self, **kwargs):
        """ Return a formated experiments query url for the calls arguments
        """
        return self.query_url("files", **kwargs)
    
    def query_experiment(self, **kwargs):
        url = self.query_url_experiments(**kwargs)
        stream = urllib2.urlopen(url, timeout=self.timeout)
        #  TODO: check stream for errors  
        return stream
    
    def query_files(self, **kwargs):
        url = self.query_url_files(**kwargs)
        stream = urllib2.urlopen(url, timeout=self.timeout)
        #  TODO: check stream for errors  
        return stream
    
    def open_file(self, accession, kind="raw", ext=None):
        """ Return a file handle to experiment data.
        Possible values for kind:
            - raw: return the raw data if available
            - fgem: return the processed data if available
            - biosamples: a png or svg design image
            - idf: investigation description
            - adf: array design description
            - mageml: MAGE-ML file
            
        Example::
        
            >>> raw_file = conn.open_file("E-TABM-1087", kind="raw")
            >>> processed_file = conn.open_file("E-TABM-1087", kind="fgem")
             
        """
        from Orange.misc.xml import parse 
        files = parse(self.query_files(accession=accession))
        files = list(files.elements("file"))
        for file in files:
            filekind = file.elements("kind").next()
            fileext = file.elements("extension").next()
            if filekind.data.strip() == kind and (fileext.data.strip() == ext or ext is None): 
                url = file.elements("url").next()
                return urllib2.urlopen(url.data.strip(), timeout=self.timeout)
    
    
def query_experiments(**kwargs):
    """ Query Array Express experiments.
    
    Example ::
    
        >>> query_experiments(species="Homo sapiens", ef="organism_part", efv="liver")
        <addinfourl at ...
        
    """
    return ArrayExpressConnection().query_experiment(**kwargs)

def query_files(**kwargs):
    """ Query Array Express files.
    
    Example ::
    
        >>> query_files(species="Mus musculus", ef="developmental_stage", efv="embryo")
        <addinfourl at ...
                        
    """
    return ArrayExpressConnection().query_files(**kwargs)
    
# TODO: List all accepted keyword values for query_* functions.

"""\
Gene Expression Atlas
---------------------
"""

class GeneExpressionAtlasConenction(object):
    """ A connection to Gene Expression Atlas database
    """
    try:
        DEFAULT_CACHE = shelve.open(os.path.join(orngEnviron.bufferDir, "GeneExpressionAtlasCache.shelve"))
    except ImportError:
        warnings.warn("Could not load persistent cache!")
        DEFAULT_CACHE = {}
    
    DEFAULT_ADDRESS = "http://www.ebi.ac.uk:80/gxa/"
    
    def __init__(self, address=None, cache=None, timeout=30):
        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.cache = cache if cache is not None else self.DEFAULT_CACHE
        self.timeout = timeout
        
    def format_query(self,):
        pass
    
    def query(self, condition, format="json", start=None, rows=None, indent=False):
        url = self.address + "api?" + condition.rest()
        if start and rows:
            url += "&start={0}&rows={1}".format(start, rows)
        url += "&format={0}".format(format)
        if indent:
            url += "&indent"
#        print url
        response = urllib2.urlopen(url)
        return response
    
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
    
GENE_FILTER_QUALIFIERS =\
    ["Is",
     "IsNot"
     ]

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
     "Unknown",
     "Xenopus laevis"
     ]
    
def ef_ontology():
    """ Return the EF (Experimental Factor) ontology
    """
    import obiOntology
#    return obiOntology.OBOOntology(urllib2.urlopen("http://efo.svn.sourceforge.net/svnroot/efo/trunk/src/efoinobo/efo.obo"))
    import orngServerFiles
    # Should this be in the OBOFoundry (Ontology) domain
    file_name = orngServerFiles.localpath_download("ArrayExpress", "efo.obo")
    return obiOntology.OBOOntology(open(filename, "rb"))


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
#        assert(self.factor in efv_ontology())
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
        
    
def query_atlas_simple(genes=None, regulated=None, organism=None, format="json"):
    """ A simple Atlas query.
    
    Example::
        
        >>> query_atlas_simple(genes=['Pou5f1', 'Dppa3'], organism="Mus musculus")
        <addinfourl at ...
        
        >>> query_atlas_simple(genes=['Pou5f1', 'Dppa3'], regulated="up", organism="Mus musculus")
        <addinfourl at ...
        
    """
    conditions = AtlasConditionList()
    conditions.append(AtlasConditionGeneProperty("Gene", "Is", genes))
    if regulated:
        conditions.append(AtlasConditionExperimentalFactor("", regulated, 1, ""))
    if organism:
        conditions.append(AtlasConditionOrganism(organism))
    connection = GeneExpressionAtlasConenction()
    results = connection.query(conditions, format=format)
    return results

"""\
TODO: can this be implemented query_atlas(organism="...", Locuslink="...", Chebi="...", up3InCompound="..." downInEFO="...")
      Need a full list of accepted factors 
"""

def query_atlas(condition, format="json", start=None, rows=None, indent=False):
    """ Query Atlas based on a `condition` (instance of AtlasCondition)
    
    Example::
        
        >>> #condition = AtlasConditionGeneProperty()
        
    """
    connection = GeneExpressionAtlasConenction()
    results = connection.query(condition, format=format, start=start,
                               rows=rows, indent=indent)
    return results


def get_atlas_summary(genes, organism):
    """ Return 3 dictionaries containing a summary of atlas information
    about three experimental factors:
    
        - Organism Part (OP) 
        - Disease State (DS)
        - Cell type (CT)
    
    Each dictionary contains query genes as keys. Values are dictionaries
    mapping factor values to a 2-tuple containig the count of up regulated
    and down regulated experiments.
    
    Example::
    
        >>> get_atlas_summary(["RUNX1"], "Homo sapiens")
        ({u'RUNX1': ...
        
    """
    genes_condition = AtlasConditionGeneProperty("Gene", "Is", genes)
    org_condition = AtlasConditionOrganism(organism)
    condition = AtlasConditionList([genes_condition, org_condition])
    result = query_atlas(condition, format="json")
    import json
    result = json.load(result)
    
    org_part = collect_ef_summary(result, "organism_part")
    disease_state = collect_ef_summary(result, "disease_state")
    cell_type = collect_ef_summary(result, "cell_type")
    
    return org_part, disease_state, cell_type
    
def collect_ef_summary(info, ef):
    """ Collect the results summary from query_atlas, result for experimental
    factor `ef`. 
    """
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
    
    return dict(summary)
    
    
if __name__ == "__main__":
    from pprint import pprint    
    pprint(get_atlas_summary(['Pou5f1', 'Dppa3'], 'Mus musculus'))
       
    pprint(get_atlas_summary(['PDLIM5', 'FGFR2' ], 'Homo sapiens'))
    
    
    conn = ArrayExpressConnection()
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS, extraglobs={"conn": conn})
    