"""
===============
obiArrayExpress
===============

A python module for accessing the ArrayExpress web services and database.

`Array Express Archive <http://www.ebi.ac.uk/arrayexpress/>`_ is a
database of gene expression experiments that you can query and download.

Example ::

    >>> # Retrieve the object representing experiment with accession E-TABM-25
    >>> experiment = ArrayExpressExperiment("E-TABM-25")
    >>> print experiment.accession
    E-TABM-25

    >>> print experiment.name
    Transcription profiling of aging in the primate brain

    >>> print experiment.species
    ['Pan troglodytes']

    >>> print experiment.files
    [{'kind': ...

    >>> # Retrieve the data matrix for experiment 'E-MEXP-2917'
    >>> experiment = ArrayExpressExperiment("E-MEXP-2917")
    >>> table = experiment.fgem_to_table()


Low level Array Express query using REST services::

    >>> from Orange.bio import obiArrayExpress
    >>> obiArrayExpress.query_experiments(accession='E-MEXP-31')
    {u'experiments': ...

    >>> obiArrayExpress.query_experiments(keywords='gliobastoma')
    {u'experiments': ...

    >>> obiArrayExpress.query_files(accession='E-MEXP-32', format="xml")
    <xml.etree.ElementTree.ElementTree object ...

.. note:: Currently querying ArrayExpress files only works with the xml format.

.. note:: See the documentation of `query_experiments` for a full set of
          parameters that these functions accept.

"""

from __future__ import absolute_import

import os
import urllib2
import re
import warnings
import shelve
import shutil
import posixpath
import json
from xml.etree.ElementTree import ElementTree
from StringIO import StringIO

from collections import defaultdict
from functools import partial
from contextlib import closing

from Orange.orng import orngServerFiles

parse_json = json.load


def parse_xml(stream):
    """ Parse an xml stream into an instance of
    `xml.etree.ElementTree.ElementTree`.
    """
    return ElementTree(file=stream)

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


class ArrayExpressConnection(object):
    """
    A connection to the Array Express.

    Used to construct a REST query and run it.

    :param address: Address of the ArrayExpress API.
    :param timeout: Timeout for the socket connection.

    .. todo::
        Implement user login (see Accessing Private Data in API docs)

    """

    DEFAULT_ADDRESS = "http://www.ebi.ac.uk/arrayexpress/{format}/v2/"
    DEFAULT_FORMAT = "json"
    DEFAULT_CACHE = orngServerFiles.localpath("ArrayExpress",
                                              "ArrayExpressCache.shelve")
    # Order of arguments in the query
    _ARGS_ORDER = ["keywords", "species", "array"]

    def __init__(self, address=None, timeout=30, cache=None,
                 username=None, password=None):
        self.address = address if address is not None else self.DEFAULT_ADDRESS
        self.timeout = timeout
        self.cache = cache if cache is not None else self.DEFAULT_CACHE
        self.username = username
        self.password = password
        
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
            # TODO check if val contains a datetime.date object, assert proper format
            return format_interval(val)
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
        stream = self._cache_urlopen(url, timeout=self.timeout)
        return stream
    
    def query_files(self, **kwargs):
        """ Return an open stream to the files query results.
        
        .. note:: This member function takes the same arguments as the module
                  level `query_files` function.
        """
        url = self.query_url_files(**kwargs)
        stream = self._cache_urlopen(url, timeout=self.timeout)
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
        stream = self.query_files(accession=accession, format="json")
        data = json.load(stream)
        try:
            files = data["files"]["experiment"]["file"]
        except KeyError:
            raise ValueError(accession)

        for file in files:
            filekind = file["kind"]
            fileext = file["extension"]
            if (filekind == kind) and (fileext == ext or ext is None):
                url = file["url"]
                return self._cache_urlopen(str(url), timeout=self.timeout)

        raise ValueError("%s does not have a file of kind: %r" %
                         (accession, kind))

    def _cache_urlopen(self, url, timeout=30):
        if self.cache is not None:
            with self.open_cache("r") as cache:
                if url in cache:
                    return StringIO(cache[url])

            stream = urllib2.urlopen(url, timeout=timeout)
            data = stream.read()
            with self.open_cache("w") as cache:
                cache[url] = data

            return StringIO(data)
        else:
            return urllib2.urlopen(url, timeout=timeout)

    def open_cache(self, flag="r"):
        if isinstance(self.cache, basestring):
            try:
                return closing(_open_shelve(self.cache, flag))
            except Exception:
                return fake_closing({})
        elif hasattr(self.cache, "close"):
            return closing(self.cache)
        elif self.cache is None:
            return fake_closing({})
        else:
            return fake_closing(self.cache)

from contextlib import contextmanager

@contextmanager
def fake_closing(obj):
    yield obj
    
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
    
    
def query_files(keywords=None, accession=None, array=None, ef=None,
               efv=None, expdesign=None, exptype=None,
               gxa=None, pmid=None, sa=None, species=None,
               expandefo=None, directsub=None, assaycount=None,
               efcount=None, samplecount=None, rawcount=None,
               fgemcount=None, miamescore=None, date=None, 
               format="json", wholewords=None, connection=None):
    """ Query Array Express files.
    
    This function accepts the same arguments as `query_experiments`.
    
    Example ::
    
        >>> query_files(species="Mus musculus", ef="developmental_stage", efv="embryo", format="xml")
        <xml.etree.ElementTree.ElementTree object ...
        
    .. todo:: why does the json interface not work.
                        
    """
    if connection is None:
        connection = ArrayExpressConnection()
    
    stream = connection.query_files(keywords=keywords, accession=accession,
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
    
    
def open_file(accession, kind="raw", ext=None, repo_dir=None):
    """ Open a file for the experiment. 
     
    Example ::
    
        >>> file = open_file("E-MTAB-369", kind="raw", repo_dir="~/.arrayexpress/")
        
    """
    raise NotImplementedError

"""\
MAGE-TAB convinence functions, classes
======================================
"""

IDF_SINGLE_TAGS = \
    ["Investigation Title",
     "Date of Experiment",
     "Public Release Date",
     "Experiment Description",
    ]

def parse_idf(file):
    """ Parse an idf.txt (Investigation Description  Format) formated file.
    Return a list of tuples where the first element is the tag (first column)
    and the second element is a list of all tag values.
    
    """
    if isinstance(file, basestring):
        file = open(file, "rb")
    data = file.read()
    lines = data.splitlines()
    lines = [line.split("\t") for line in lines if line and not line.startswith("#")]
    parsed = [(line[0], line[1:]) for line in lines]
    return parsed
    
def parse_sdrf(file):
    """ Parse an sdfr formated file. Return a tuple with the first element
    a list of header values and the second element is a list of all rows
    (a row is a list of string values).
    
    """
    if isinstance(file, basestring):
        file = open(file, "rb")
    data = file.read()
    lines = data.splitlines()
    lines = [line.split("\t") for line in lines if line.strip() and not line.startswith("#")]
    header = [h for h in lines[0] if h]
    rows = [line[:len(header)] for line in lines[1:]]
    assert(all([len(r) == len(header) for r in rows]))
    return header, rows
    
def parse_adf(file):
    pass

def parse_data_matrix(file):
    """ Parse the MAGE-TAB processed data matrix. Return a tuple where the
    elements are:
        - a (header REF, header values) tuple (e.g. ("Hybridization REF", ["Stage1", "Stage2", ...]) )
        - a list of quantitation type for header values (e.g. ["log2 ratio", "log2 ratio", ...])
        - a (row REF, row names list) tuple ("e.g. ("Reporter REF", ["Probe 1", "Probe 2", ...]) )
        - a list of list matrix with values (as strings)
        
    """
    if isinstance(file, basestring):
        file = open(file, "rb")
    data = file.read()
    lines = data.splitlines()
    lines = [line.split("\t") for line in lines if line.strip()]
    header = lines[0]
    header_ref, header = header[0], header[1:]
    line2 = lines[1]
    row_ref, quant_type = line2[0], line2[1:]
    row_names, rows = [], []
    for line in lines[2:]:
        r_name, row = line[0], line[1:]
        row_names.append(r_name)
        rows.append(row)
        
    return ((header_ref, header),
            quant_type,
            (row_ref, row_names),
            rows) 
    
class InvestigationDesign(dict):
    """ Investigation design (contains the contents of the .idf).
    
    Example ::
    
        >>> idf = InvestigationDesign("foobar.idf")
        >>> print idf.investigation_title
        foo investigation
        >>> print idf.experimental_design
        ['fubar', 'snafu']
        >>> print idf.sdrf_file
        ['foobar.sdrf']
        
    """
    def __init__(self, idf_file=None):
        idf = parse_idf(idf_file)
        self._idf = idf
        self.update(dict(idf))
        for tag, values in idf:
            if tag in IDF_SINGLE_TAGS:
                values = values[0] if values else None
            ttag = self.transform_tag(tag)
            setattr(self, ttag, values)
        
    def transform_tag(self, tag):
        """ Transform the tag into a proper python attribute name by
        replacing all spaces and special characters (e.g '[', ']' into
        underscores).
        
        """
        toreplace = [" ", "-", "[", "]"]
        for s in toreplace:
            tag = tag.replace(s, "_")
        return tag.lower()
            
    def __getitem__(self, tag):
        """ Return the tag values
        
        Example ::
        
            >>> idf["Investigation Title"]
            'foo investigation'
            
        """
        try:
            return self._idf_dict[tag]
        except KeyError:
            pass
        
        ttag = self.transform_tag(tag)
        if hasattr(self, ttag):
            return getattr(self, ttag)
        else:
            raise KeyError(tag)
        
class SampleDataRelationship(object):
    """ Sample-Data Relationship (contains the contents of the .sdrf file).
    
    Example ::
    
        >>> sdr = SampleDataRelationship("foobar.sdrf")
        >>> sdr.source_name()
        ['foo', ...
        
        >>> sdr.sample_name()
        ['sampled foo', ...
        
        >>> sdr.extract_protocol_ref()
        ['bar the foo', ...
        
        >>> sdr.derived_array_data_matrix_file()
        ['foobar.data.txt', ...
        
        >>>
         
    """
    
    # Nodes an edges
    NODES_EDGES = ["Source Name", "Sample Name", "Extract Name",
                   "Labeled Extract Name", "Hybridization Name",
                   "Assay Name", "Scan Name", "Normalization Name",
                   "Array Data File", "Derived Array Data File",
                   "Array Data Matrix File", "Derived Array Data Matrix File",
                   "Image File", "Protocol REF"]
    
    # Attributes for nodes and edges
    NODE_EDGE_ATTRIBUTES = \
    {"Source Name": ["Characteristics", "Provider", "Material Type", "Description", "Comment"],
     "Sample Name": ["Characteristics", "Material Type", "Description", "Comment"],
     "Extract Name":["Characteristics", "Material Type", "Description", "Comment"],
     "Labeled Extract Name": ["Characteristics", "Material Type", "Description", "Label", "Comment"],
     "Hybridization Name": ["Array Design File", "Array Design REF", "Comment"],
     "Assay Name": ["Technology Type", "Comment"],
     "Scan Name": ["Comment"],
     "Normalization Name": ["Comment"],
     "Array Data File": ["Comment"],
     "Derived Array Data File": ["Comment"],
     "Array Data Matrix File": ["Comment"],
     "Derived Array Data Matrix File": ["Comment"],
     "Image File": ["Comment"],
     "Protocol REF": ["Term Source REF", "Parameter", "Performer", "Date", "Comment"]
     }
    
    # Attributes
    ATTRIBUTE_COLUMNS = \
    {"Characteristics []": ["Unit", "Term Source REF"],
     "Provider": ["Comment"],
     "Material Type": ["Term Source REF"],
     "Label": ["Term Source REF"],
     "Array Design File": ["Comment"],
     "Array Design REF": ["Term Source REF", "Comment"],    
     "Technology Type": ["Term Source REF"],
     "Factor Value [] ()": ["Unit", "Term Source REF"],
     "Performer": ["Comment"],
     "Date": [],
     "Parameter Value []": ["Unit", "Comment"],
     "Unit []": ["Term Source REF"],
     "Description": [],
     "Term Source REF": ["Term Accession Number"],
     "Term Accession Number": [],
     "Comment []": []
     }
    def __init__(self, sdrf_file=None):
        header, rows = parse_sdrf(sdrf_file)
        self.header = header
        self.rows = rows 
        
    def transform_tag(self, tag):
        """ Transform the tag into a proper python attribute name by
        replacing all spaces and special characters (e.g '[', ']' into
        underscores).
        
        """
        toreplace = [" ", "-", "[", "]"]
        for s in toreplace:
            tag = tag.replace(s, "_")
        return tag.lower()
    
    def _subsection(self, name):
        """ Return the named subsection (name must be from the
        NODES_EDGES list).
        
        """
        idx = self.NODES_EDGES.index(name)
        start = self.header.index(name)
        end = -1
        for name in self.NODES_EDGES[idx + 1:]:
            if name in self.header[start + 1:]:
                end = self.header.index(name, start + 1)
                break
        return self.header[start:end], [r[start:end] for r in self.rows]
    
    def _column(self, name):
        """ Return the named column.
        """
        if isinstance(name, basestring):
            index = self.header.index(name)
        else:
            index = name
        return [r[index] for r in self.rows]
        
    def source(self):
        """ Return the Source subsection
        """
        return self._subsection("Source Name")
        
    def source_name(self):
        """ Return the Source Name subsection
        """
        return self._column("Source Name")
        
    def sample(self):
        """ Return the Sample subsection
        """
        return self._subsection("Sample Name")
        
    def sample_name(self):
        """ Return the Sample Name subsection
        """
        return self._column("Sample Name")
        
    def extract(self):
        """ Return the Extract subsection
        """
        return self._subsection("Extract Name")
        
    def extract_name(self):
        """ Return the Extract Name subsection
        """
        return self._column("Extract Name")
        
    def labeled_extract(self):
        """ Return the Labeled Extract subsection
        """
        return self._subsection("Labeled Extract Name")
        
    def labeled_extract_name(self):
        """ Return the Labeled Extract Name subsection
        """
        return self._column("Labeled Extract Name")
        
    def hybridization(self):
        """ Return the Hybridization subsection.
        """
        return self._subsection("Hybridization Name")
        
    def hybridization_name(self):
        """ Return the Hibridization Name subsection.
        """
        return self._column("Hybridization Name")
        
    def assay(self):
        """ Return the Assay subsection
        """
        return self._subsection("Assay Name")
        
    def assay_name(self):
        """ Return the Assay Name subsection
        """
        return self._column("Assay Name")
        
    def scan(self):
        """ Return the Scan subsection
        """
        return self._subsection("Scan Name")
        
    def scan_name(self):
        """ Return the Scan name subsection
        """
        return self._column("Scan Name")
        
    def normalization(self):
        """ Return the Normalization subsection.
        """
        return self._subsection("Normalization Name")
        
    def normalization_name(self):
        """ Return the Normalization Name subsection.
        """
        return self._column("Normalization Name")
         
    def array_data(self):
        """ Return the Array Data subsection
        """
        return self._subsection("Array Data File")
    
    def array_data_file(self):
        """ Return the Array Data File subsection
        """
        return self._column("Array Data File")
        
    def derived_array_data(self):
        """ Return the Derived Array Data subsection
        """
        return self._subsection("Derived Array Data File")
    
    def derived_array_data_file(self):
        """ Return the Derived Array Data File subsection
        """
        return self._column("Derived Array Data File")
        
    def array_data_matrix(self):
        """ Return the Array Data Matrix subsection.
        """
        return self._subsection("Array Data Matrix File")
    
    def array_data_matrix_file(self):
        """ Return the Array Data Matrix File subsection.
        """
        return self._column("Array Data Matrix File")
        
    def derived_array_data_matrix(self):
        """ Return the Derived Array Data Matrix subsection.
        """
        return self._subsection("Derived Array Data Matrix File")
    
    def derived_array_data_matrix_file(self):
        """ Return the Derived Array Data Matrix File subsection.
        """
        return self._column("Derived Array Data Matrix File")
        
    def image(self):
        """ Return the Image subsection
        """
        return self._subsection("Image File")
    
    def image_file(self):
        """ Return the Image File subsection.
        """
        return self._column("Image File")
        
class ArrayDesign(object):
    """ Arary design (contains the contents of the .adf file).
    """
    def __init__(self, adf_file=None):
        adf = parse_adf(adf_file)
        self._adf = adf
    
def _is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False
    
def _is_continuous(items, check_count=100):
    """ Are the strings in items continous numbers. 
    """
    count = 0
    i = 0
    for i, item in enumerate(items):
        if _is_float(item):
            count += 1
        if i >= check_count:
            break
    return count >= i * 0.5
    
def processed_matrix_to_orange(matrix_file, sdrf=None):
    """ Load a single processed matrix file in to an Orange.data.Table
    instance. 
    """
    import numpy
    import Orange
    
    if isinstance(matrix_file, basestring):
        matrix_file = open(matrix_file, "rb")

    header, quant_type, rows, matrix = parse_data_matrix(matrix_file)
    header_ref, header = header
    row_ref, rows = rows
    
    matrix = numpy.array(matrix, dtype=object)
    
    
    features = []
    is_float = numpy.frompyfunc(_is_float, 1, 1) # an numpy ufunc
         
    for header_name, quant, column in zip(header, quant_type, matrix.T):
        if _is_continuous(column):
            feature = Orange.feature.Continuous(header_name)
            column[numpy.where(1 - is_float(column))] = "?" # relace all non parsable floats with '?'
        else:
            values = set(column)
            feature = Orange.feature.Discrete(header_name, values=sorted(values))
        feature.attributes["quantitation type"] = quant
        features.append(feature)
        
    row_ref_feature = Orange.feature.String(row_ref)
    domain = Orange.data.Domain(features, None)
    domain.addmeta(Orange.feature.Descriptor.new_meta_id(), row_ref_feature)
    
    table = Orange.data.Table(domain, [list(row) for row in matrix])
    table.setattr("header_ref", header_ref)
    # Add row identifiers
    for instance, row in zip(table, rows):
        instance[row_ref_feature] = row

    if sdrf is not None:
        pattern = re.compile(r"((Characteristics)|(Factor Value)|(Parameter Value)) \[(?P<name>.*?)\].*")
        # Row name in sdrf
        row_name = header_ref[:header_ref.find(" REF")] + " Name"
        # feature names as listed in sdrf
        feature_names = sdrf._column(row_name)
        annotations = defaultdict(partial(defaultdict, set))
        for i, header in enumerate(sdrf.header):
            match = pattern.match(header)
            if match:
                name = match.group("name")
                for val, feature_name in zip(sdrf._column(i), feature_names):
                    annotations[feature_name][name].add(val)

        def to_str(values):
            if len(values) > 1:
                return str(list(values))
            else:
                return str(list(values)[0])

        for feature in table.domain.features:
            feature.attributes.update([(key, to_str(value)) for key, value in \
                                       annotations[feature.name].items()])
    return table


def mage_tab_to_orange(idf_filename):
    """Convert an MAGE-TAB annotated experiment into an Orange.data.Table
    instance (assumes all the associated MAGE-TAB files are in the same
    directory.

    """
    dirname = os.path.dirname(idf_filename)
    idf = InvestigationDesign(idf_filename)
    
    sdrf_filename = os.path.join(dirname, idf.sdrf_file[0])
    sdrf = SampleDataRelationship(sdrf_filename)
    
    data_matrices = set(sdrf.derived_array_data_matrix_file())
    data_matrices = [name for name in data_matrices if name.strip()]
    
    tables = []
    for filename in data_matrices:
        matrix_file = os.path.join(dirname, filename)
        table = processed_matrix_to_orange(matrix_file, sdrf)
        tables.append(table)
    table = hstack_tables(tables)

    return table


def hstack_tables(tables):
    """ Stack the tables horizontaly.
    """
    import Orange
    max_len = max([len(table) for table in tables])
    stacked_features = []
    stacked_values = [[] for i in range(max_len)]
    stacked_meta_features = []
    stacked_meta_values = [{} for i in range(max_len)]
    
    for table in tables:
        stacked_features.extend(table.domain.variables)
        stacked_meta_features.extend(table.domain.getmetas().items())
        
        for i, instance in enumerate(table):
            stacked_values[i].extend(list(instance))
            stacked_meta_values[i].update(instance.getmetas())
            
        # Fill extra lines with unknowns
        for i in range(len(table), max_len):
            stacked_values[i].extend(["?"] * len(table.domain.variables))
        
    domain = Orange.data.Domain(stacked_features, tables[-1].domain.class_var)
    domain.addmetas(dict(set(stacked_meta_features)))
    table = Orange.data.Table(domain, stacked_values)
    
    # Add meta attributes
    for instance, metas in zip(table, stacked_meta_values):
        for m, val in metas.iteritems():
            instance[m] = val
            
    return table
    
def _dictify(element):
    """ 'Dictify' and xml.etree.Element.Element instance by taking all
    subnode tags as keys and the corresponding text values as items i.e. turn
    `<node><bla>foo</bla><a>b</b></node>` into a {"bla":"foo", "a":b"}
        
    """
    if element is None:
        element = []
    dict = {}
    strip = lambda s: s.strip() if s else s
    for node in element:
        dict[node.tag] = strip(getattr(node, "text", None))
    return dict
    
class SearchableList(list):
    """ A list with an advanced `search` method
    """
    def search(self, **kwargs):
        """ Search the list for elements with members that correspond the
        supplied keyword arguments.
        
            >>> foo.bar = "foo"
            >>> list = SearchableList([foo, bar])
            >>> list.search(bar="foo") # Search for objects which have a member named "bar" and that member equals "foo"
            [<__main__.foo object ...

        """
        ret = []
        for obj in self:
            for member, value in kwargs.items():
                if hasattr(obj, member) and getattr(obj, member) == value:
                    ret.append(obj)
        return ret
    
class ArrayExpressExperiment(object):
    """ An convinience class representing an Array Express Experiment.
    
    Example ::
    
        >>> ae = ArrayExpressExperiment("E-MEXP-2917")
        >>> print ae.name
        Characterization of Data Variation in Gene Expression Profiling of Human Peripheral Blood Samples
        
        >>> for file in ae.files:
        ...     print file["name"], file["url"]
        E-MEXP-2917.sdrf.txt http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-2917/E-MEXP-2917.sdrf.txt
        ...
        
        >>> table = ae.fgem_to_table() # Retieve the experiment data table 
            
    """
    
    def __init__(self, accession, connection=None):
        self.accession = accession
        self.connection = connection
        self._etree = tree = query_experiments(accession=accession, connection=self.connection, format="xml")
        experiments = tree.findall("experiment")
        # find the exact match (more then one experiment can be listed in the query result)
        experiments = [e for e in experiments if e.find("accession").text.strip() == accession]
        self._experiment = experiment = experiments[0]
        
        self.species = [e.text for e in experiment.findall("species")]
        bool_values = {"true": True, "false": False}
        self.rawdatafiles = bool_values[experiment.find("rawdatafiles").get("available","false")]
        self.fgemdatafiles = bool_values[experiment.find("processeddatafiles").get("available", "false")]
        
        self.sampleattributes = []
        for sa in experiment.findall("sampleattribute"):
            category = sa.find("category").text.strip()
            values = [val.text for val in sa.findall("value")]
            self.sampleattributes.append((category, values))
            
        self.experimentalfactors = []
        for ef in experiment.findall("experimentalfactor"):
            name = ef.find("name").text.strip()
            values = [val.text.strip() for val in ef.findall("values")]
            self.experimentalfactors.append((name, values))
            
        self.miamescores = _dictify(experiment.find("miamescores"))
            
        self.id = experiment.find("id").text
        self.secondaryaccession = getattr(experiment.find("secondaryaccession"), "text", None)
        self.name = experiment.find("name").text
        self.experimenttype = experiment.find("experimenttype").text.strip()
        self.releasedate = experiment.find("releasedate").text
        self.lastupdatedate = getattr(experiment.find("lastupdatedate"), "text", None)
        self.samples = int(experiment.find("samples").text)
        self.assays = int(experiment.find("assays").text)
        
        self.arraydesign = [_dictify(e) for e in experiment.findall("arraydesign")]
            
        self.bioassaydatagroups = [_dictify(group) for group in experiment.findall("bioassaydatagroup")]
        self.bibliography = [_dictify(e) for e in experiment.findall("bibliography")]
        self.provider = [_dictify(e) for e in experiment.findall("provider")]
        
        self.experimentdesign = []
        for expd in experiment.findall("experimentdesign"):
            self.experimentdesign.append(expd.text)
            
        self.description = [_dictify(e) for e in experiment.findall("description")]
        
        tree = query_files(accession=self.accession, format="xml", connection=self.connection)
        experiments = tree.findall("experiment")
        experiments = [e for e in experiments if e.find("accession").text.strip() == accession]
        experiment = experiments[0]
        files = experiment.findall("file")
        self.files = [_dictify(file) for file in files]
        
    def _download_processed(self):
        """ Download the processed matrix file, and associated MAGE-TAB files (idf, sdrf, adf)
        
        .. todo:: What about the raw data files (we need converters for other file types) 
        """
        assert(self.fgemdatafiles)
        exp_files = [(f["kind"], f) for f in self.files if f.get("kind") in ["idf", "sdrf"] and f.get("extension") == "txt"]
        exp_files += [(f["kind"], f) for f in self.files if f.get("kind") == "fgem"]
        array_files = [(f["kind"], f) for f in self.files if f.get("kind") == "adf" and f.get("extension") == "txt"]
        assert(len(files) == 3)
        
        for type, file in files.iteritems():
            url = file["url"].strip()
            rest, basename = os.path.split(url)
            _, dirname = os.path.split(rest)
            
            repo_dir = orngServerFiles.localpath("ArrayExpress", dirname)
            try:
                os.makedirs(repo_dir)
            except OSError:
                pass
            local_filename = os.path.join(repo_dir, basename)
            stream = urllib2.urlopen(url)
            shutil.copyfileobj(stream, open(local_filename, "wb"))
            
            if file["extension"] == "zip":
                import zipfile
                zfile = zlib.ZipFile(local_filename)
                zfile.extractall(repo_dir)
            elif file["extension"] == "gz":
                import gzip
                gzfile = gzip.open(local_filename)
                gzfile.extractall(repo_dir)
            elif file["extension"] in ["tgz", "tar"]:
                import tarfile
                tfile = tarfile.TarFile(local_filename)
                tfile.extractall(repo_dir)
            elif file["extension"] == "txt":
                pass
            else:
                raise ValueError("Unknown extension ('{0}').".format(basename))
            
    def _download_file(self, url, extract=True):
        """ Download the `file` from the ArrayExpress saving it to a local
        repository directory.
         
        """
        rest, basename = posixpath.split(url)
        dirname = posixpath.basename(rest)
        repo_dir = orngServerFiles.localpath("ArrayExpress", dirname)
        try:
            os.makedirs(repo_dir)
        except OSError:
            pass
        stream = urllib2.urlopen(url)
        local_filename = os.path.join(repo_dir, basename)
        shutil.copyfileobj(stream, open(local_filename, "wb"))
        
        if extract:
            _, extension = os.path.splitext(local_filename)
            if extension == ".zip":
                import zipfile
                zfile = zipfile.ZipFile(local_filename)
                zfile.extractall(repo_dir)
            elif extension == ".gz":
                import gzip
                gzfile = gzip.open(local_filename)
                gzfile.extractall(repo_dir)
            elif extension in [".tgz"]:
                import tarfile
                tfile = tarfile.TarFile(local_filename)
                tfile.extractall(repo_dir)
            elif extension == ".txt":
                pass
            else:
                raise ValueError("Unknown extension ('{0}').".format(basename))
            
    def _is_local(self, url):
        """ Is the `url` stored in the local repository.
        """
        return os.path.exists(self._local_filepath(url))
    
    def _local_filepath(self, url):
        """ Return the local file path for url.
        """
        rest, basename = posixpath.split(url)
        dirname = posixpath.basename(rest)
        return orngServerFiles.localpath("ArrayExpress", os.path.join(dirname, basename))
    
    def _open(self, url):
        """ Return an open file like handle to url (ArrayExpress file).
        The file is cached in the local repository for future access.
        
        """
        if not self._is_local(url):
            self._download_file(url, extract=True)
        file = self._local_filepath(url)
        return open(file, "rb")
    
    def _search_files(self, kind=None, extension=None):
        """ Search files by `kind` and `extension`.
        """
        res = []
        for file in self.files:
            kind_match = kind == file.get("kind") or kind is None
            extension_match = extension == file.get("extension") or extension is None
            
            if kind_match and extension_match:
                res.append(file)
        return res
        
    def array_design(self):
        """ Return a list of `ArrayDesign` instances used in this experiment.
        """
        files = [f for f in self.files if f.get("kind") == "adf" and \
                 f.get("extension") == "txt"]
        
        array_design = []
        for file in files:
            url = file.get("url")
            if not self._is_local(url):
                self._download_file(url)
            array_design.append(ArrayDesign(self._open(url)))
        return array_design
        
    def investigation_design(self):
        """ Return an `InvestigationDesgin` instance for this experiment
        """
        files = [f for f in self.files if f.get("kind") == "idf" and \
                 f.get("extension") == "txt"]
        if not files:
            raise ValueError("The experiment '{0}' does not have an investigation design file".format(self.accession))
        file = files[0]
        return InvestigationDesign(self._open(file.get("url")))
        
        
    def sample_data_relationship(self):
        """ Return an `SampleDataRelationship` instance describing this experiment.
        """
        files = [f for f in self.files if f.get("kind") == "sdrf" and \
                 f.get("extension") == "txt"]
        if not files:
            raise ValueError("The experiment '{0}' does not have an sample and data relationship file".format(self.accession))
        file = files[0]
        return SampleDataRelationship(self._open(file.get("url")))
        
    def fgem_to_table(self):
        """ Retrieve the processed matrix from the Array Express ftp
        server and convert it to an :class:`Orange.data.Table` instance.
         
        """
        assert(self.fgemdatafiles)
        repo_dir = orngServerFiles.localpath("ArrayExpress", self.accession)
        # Find the file listing the data matrix files (should be in sdrf but sometimes it is in 2column file only, why?)
        sdrf = self._search_files("sdrf", "txt")
        if sdrf:
            sdrf = SampleDataRelationship(self._open(sdrf[0].get("url")))
            if "Derived Array Data Matrix File" not in sdrf.header:
                twocol = self._search_files("twocolumn", "txt")
                if twocol:
                    sdrf = SampleDataRelationship(self._open(twocol[0].get("url")))
        matrix_file = self._search_files("fgem")[0]
        self._open(matrix_file.get("url")) 
        matrix_files = sorted(set(sdrf.derived_array_data_matrix_file()))
        
        idf_file = self._search_files("idf", "txt")[0]
        self._open(idf_file.get("url")) # To download if not cached
        return mage_tab_to_orange(os.path.join(repo_dir, idf_file.get("name")))
        
    
"""\
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
    DEFAULT_ADDRESS = "http://www.ebi.ac.uk/gxa/api/deprecated"

    DEFAULT_CACHE = orngServerFiles.localpath(
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
                return fake_closing({})
        else:
            return fake_closing(self.cache)


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
##     "Unknown",
#     "Xenopus laevis"
#     }
    
def ef_ontology():
    """ Return the `EF <http://www.ebi.ac.uk/efo/>`_ (Experimental Factor) ontology
    """
    from . import obiOntology
#    return obiOntology.OBOOntology(urllib2.urlopen("http://efo.svn.sourceforge.net/svnroot/efo/trunk/src/efoinobo/efo.obo"))
    from Orange.orng import orngServerFiles
    # Should this be in the OBOFoundry (Ontology) domain
    try:
        file = open(orngServerFiles.localpath_download("ArrayExpress", "efo.obo"), "rb")
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

"""\
.. todo:: can this be implemented query_atlas(organism="...", Locuslink="...", Chebi="...", up3InCompound="..." downInEFO="...")
      Need a full list of accepted factors 
"""

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


def test():    
    conn = ArrayExpressConnection()
    import doctest
    foo = type("foo", (object,), {})()
    bar = type("bar", (object,), {})()
    
    doctest.testmod(optionflags=doctest.ELLIPSIS,
                    extraglobs={"conn": conn, "foo": foo,"bar": bar})
    
    
if __name__ == "__main__":
    test()
    
