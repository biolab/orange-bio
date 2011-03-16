"""
obiArrayExpress
===============

A python module for accessing the ArrayExpress and GeneExpressionAtlas
web services.


Array Express
-------------

`Array Express Archive <http://www.ebi.ac.uk/arrayexpress/>`_ is a database of gene expression experiments that you
can query and download.

Example of an Array Express query ::

    >>> import obiArrayExpress
    >>> obiArrayExpress.query_experiments(accession='E-MEXP-31')
    {u'experiments': ...
    
    >>> obiArrayExpress.query_files(accession='E-MEXP-32', format="xml")
    XMLNode(...
   
.. note:: Currently querying ArrayExpress files only works with the xml format.

.. note:: See the documentation of `query_experiments` for a full set of
          parameters that these functions accept.

"""

import os, sys
import urllib2

import orngEnviron
import warnings
import posixpath
import shelve
import json
from collections import defaultdict

parse_json = json.load

def parse_xml(stream):
    #TODO: parse stream, return the same structure as parse_json            
    from Orange.misc.xml import parse
    tree = parse(stream)
    return tree


# All searchable fields of ArrayExpress (see query_experiments docstring
# for a description of the fields)
ARRAYEXPRESS_FIELDS = \
    ["keywords",
     "accession",
     "array",
     "ef",
     "efv",
     "expdesign",
     "exptype",
     "gxa",
     "pmid",
     "sa",
     "species",
     "expandefo",
     "directsub",
     "assaycount",
     "efcount",
     "samplecount",
     "sacount",
     "rawcount",
     "fgemcount",
     "miamescore",
     "date",
     "wholewords",
    ]

class ArrayExpressConnection(object):
    """ A connection to the ArrayExpress. Used for query construction,
    and user login. 
    
    .. todo:: Implement user login.
    """
    
    DEFAULT_ADDRESS = "http://www.ebi.ac.uk/arrayexpress/{format}/v2/"
    DEFAULT_FORMAT = "json"
    
    # Order of arguments in the query
    _ARGS_ORDER = ["keywords", "species", "array"]
    
    def __init__(self, address=None, timeout=30,
                 username=None, password=None):
        """ Initialize the connection object.
        
        :param address: Address of the ArrayExpress API
        :param timeout: Timeout for the socket connection
        
        .. todo:: Implement user login (see Accessing Private Data in API docs)
        
        """
        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.timeout = timeout
        
    def format_query(self, **kwargs):
        """ Format the query arguments.
        
        Example ::
        
            >>> conn.format_query(gxa=True, efcount=(1, 5))
            'efcount=[1 TO 5]&gxa=true'
            
        """
        # Formaters:
        def format_default(val):
            if isinstance(val, basestring):
                return val
            else:
                return "+".join(val)
        def format_species(val):
            return '"%s"' % val.lower()
        def format_gxa(val):
            if val:
                return "true"
            else:
                raise ValueError("gxa={0}".format(val))
        def format_expandefo(val):
            if val:
                return "on"
            else:
                raise ValueError("expandefo={0}".format(val))
        def format_true_false(val):
            return "true" if val else "false"
        def format_interval(val):
            if isinstance(val, tuple):
                return "[{0} TO {1}]".format(*val)
            else:
                raise ValueError("Must be an interval argument (min, max)!")
        def format_date(val):
            return val
        def format_wholewords(val):
            if val:
                return "on"
            else:
                raise ValueError("wholewords={0}".format(val))
        
        formaters = {"species": format_species,
                     "gxa": format_gxa,
                     "expandefo": format_expandefo,
                     "directsub": format_true_false,
                     "assaycount": format_interval,
                     "efcount": format_interval,
                     "samplecount": format_interval,
                     "sacount": format_interval,
                     "rawcount": format_interval,
                     "fgemcount": format_interval,
                     "miamescore": format_interval,
                     "date": format_date,
                     "wholewords": format_wholewords,
                     }
        parts = []
        arg_items = kwargs.items()
        ordered = sorted(arg_items, key=lambda arg: self._ARGS_ORDER.index(arg[0]) \
                         if arg[0] in self._ARGS_ORDER else 100)
        
        for key, value in kwargs.iteritems():
            if key == "format":
                continue # format is handled in query_url
            if key not in ARRAYEXPRESS_FIELDS:
                raise ValueError("Invalid argument name: '{0}'".format(key))
            if value is not None and value != []:
                fmt = formaters.get(key, format_default)
                value = fmt(value)
                parts.append("{0}={1}".format(key, value))
                 
        return "&".join(parts)
        
    def query_url(self, what="experiments", **kwargs):
        """ Return a formated query URL for the query arguments
        
        Example ::
            >>> conn.query_url(accession="E-MEXP-31")
            'http://www.ebi.ac.uk/arrayexpress/json/v2/experiments?accession=E-MEXP-31'
            
        """
        query = self.format_query(**kwargs)
        url = posixpath.join(self.address, what)
        url = url.format(format=kwargs.get("format", self.DEFAULT_FORMAT))
        url = url + ("?" + query if query else "")
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
        """ Return an open stream to the experiments query results.
        
        .. note:: This member function takes the same arguments as the module
                  level `query_experiemnts` function.
          
        """
        url = self.query_url_experiments(**kwargs)
        stream = urllib2.urlopen(url, timeout=self.timeout)
        return stream
    
    def query_files(self, **kwargs):
        """ Return an open stream to the files query results.
        
        .. note:: This member function takes the same arguments as the module
                  level `query_files` function.
        """
        url = self.query_url_files(**kwargs)
        stream = urllib2.urlopen(url, timeout=self.timeout)
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
            
        Example ::
        
            >>> raw_file = conn.open_file("E-TABM-1087", kind="raw")
            >>> processed_file = conn.open_file("E-TABM-1087", kind="fgem")
             
        """
        from Orange.misc.xml import parse 
        files = parse(self.query_files(accession=accession), format="xml")
        files = list(files.elements("file"))
        for file in files:
            filekind = file.elements("kind").next()
            fileext = file.elements("extension").next()
            if filekind.data.strip() == kind and (fileext.data.strip() == ext or ext is None): 
                url = file.elements("url").next()
                return urllib2.urlopen(url.data.strip(), timeout=self.timeout)
    
    
def query_experiments(keywords=None, accession=None, array=None, ef=None,
                      efv=None, expdesign=None, exptype=None,
                      gxa=None, pmid=None, sa=None, species=None,
                      expandefo=None, directsub=None, assaycount=None,
                      efcount=None, samplecount=None, rawcount=None,
                      fgemcount=None, miamescore=None, date=None, 
                      format="json", wholewords=None, connection=None):
    """ Query Array Express experiments.
    
    :param keywords: A list of keywords to search (e.g. ['gliobastoma']
    :param accession: Search by experiment accession (e.g. 'E-MEXP-31')
    :param array: Search by array design name or accession (e.g. 'A-AFFY-33')
    :param ef: Experimental factor (names of main variables of experiments)
    :param efv: Experimental factor value (Has EFO expansion)
    :param expdesign: Experiment design type. (e.g. ["dose", "response"])
    :param exptype: Experiment type (e.g. 'RNA-Seq', has EFO expansion)
    :param gxa: If True limit the results to experiments from the Gene
        Expreission Atlas only.
    :param pmid: Search by PubMed identifier
    :param sa: Sample attribute values (e.g. 'fibroblast', has EFO expansion)
    :param species: Search by species (e.g. 'Homo sapiens', has EFO expansion)
    
    :param expandefo: If True expand the search terms with all its child terms
        in the Experimental Factor Ontology (EFO_) (e.g. keywords="cancer"
        will be expanded to include for synonyms and sub types of cancer).
    :param directsub: If True return only experiments submited directly to
        Array Express else if False return only experiments imported from GEO
        database (default None, return both)
    :param assaycount: A two tuple (min, max) for filter on the number of
        assays (e.g. (1, 5) will return only experiments with at least one
        and no more then 5 assays).
    :param efcount: Filter on the number of experimental factors (e.g. (1, 5))
    :param sacount: Filter on the number of sample attribute categories
    :param rawcount: Filter on the number or raw files
    :param fgemcount: Filter on the number of final gene expression matrix
        (processed data) files
    :param miamescore: Filter on the MIAME complience score (max 5)
    :param date: Filter by release date
    
    Example ::
    
        >>> query_experiments(species="Homo sapiens", ef="organism_part", efv="liver")
        {u'experiments': ...
        
    .. _EFO: http://www.ebi.ac.uk/efo/
    
    """
    if connection is None:
        connection = ArrayExpressConnection()
        
    stream = connection.query_experiment(keywords=keywords, accession=accession,
                array=array, ef=ef, efv=efv, expdesign=expdesign, exptype=exptype,
                gxa=gxa, pmid=pmid, sa=sa, species=species, expandefo=expandefo,
                directsub=directsub, assaycount=assaycount, efcount=efcount,
                samplecount=samplecount, rawcount=rawcount, fgemcount=fgemcount,
                miamescore=miamescore, date=date,  format=format,
                wholewords=wholewords)
    
    if format == "json":
        return parse_json(stream)
    else:
        return parse_xml(stream)

def query_files(**kwargs):
    """ Query Array Express files.
    
    This function accepts the same arguments as `query_experiments`.
    
    Example ::
    
        >>> query_files(species="Mus musculus", ef="developmental_stage", efv="embryo", format="xml")
        XMLNode(...
                        
    """
    connection = kwargs.get("connection", None)
    if connection is None:
        connection = ArrayExpressConnection()
    
    stream = connection.query_files(**kwargs)
    if kwargs.get("format", "json") == "json":
        return parse_json(stream)
    else:
        return parse_xml(stream)
    
__doc__ += """\
Gene Expression Atlas
---------------------

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

class GeneExpressionAtlasConenction(object):
    """ A connection to Gene Expression Atlas database.
    """
    DEFAULT_ADDRESS = "http://www.ebi.ac.uk:80/gxa/"
    
    def __init__(self, address=None, timeout=30):
        """ Initialize the conenction.
        
        :param address: Address of the server.
        :timeout: Socket timeout.
        
        """
        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.timeout = timeout
    
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
    """ Return the `EF <http://www.ebi.ac.uk/efo/>`_ (Experimental Factor) ontology
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

"""\
.. todo:: can this be implemented query_atlas(organism="...", Locuslink="...", Chebi="...", up3InCompound="..." downInEFO="...")
      Need a full list of accepted factors 
"""

def query_atlas(condition, format="json", start=None, rows=None, indent=False):
    """ Query Atlas based on a `condition` (instance of AtlasCondition)
    
    Example ::
        
        >>> condition1 = AtlasConditionGeneProperty("Goterm", "Is", "p53 binding")
        >>> condition2 = AtlasConditionExperimentalFactor("Organism_part", "up", 3, "heart")
        >>> condition = AtlasConditionList([condition1, condition2])
        >>> query_atlas(condition)
        {u'...
        
    """
    connection = GeneExpressionAtlasConenction()
    results = connection.query(condition, format=format, start=start,
                               rows=rows, indent=indent)
    if format == "json":
        return parse_json(results)
    else:
        return parse_xml(results)


def get_atlas_summary(genes, organism):
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
    genes_condition = AtlasConditionGeneProperty("Gene", "Is", genes)
    org_condition = AtlasConditionOrganism(organism)
    condition = AtlasConditionList([genes_condition, org_condition])
    result = query_atlas(condition, format="json")
    
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
    
    
def test():
    from pprint import pprint    
    pprint(get_atlas_summary(['Pou5f1', 'Dppa3'], 'Mus musculus'))
       
    pprint(get_atlas_summary(['PDLIM5', 'FGFR2' ], 'Homo sapiens'))
    
    
    conn = ArrayExpressConnection()
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS, extraglobs={"conn": conn})
    
if __name__ == "__main__":
    test()
    