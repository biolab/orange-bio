from __future__ import absolute_import

import os
import urllib2
import re
import shelve
import shutil
import posixpath
import json
from xml.etree.ElementTree import ElementTree
from StringIO import StringIO

from collections import defaultdict
from functools import partial
from contextlib import closing
from contextlib import contextmanager

from Orange.utils import serverfiles

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
    Constructs and runs REST query on ArrayExpress.

    :param address: Address of the ArrayExpress API.
    :param timeout: Timeout for the connection.

    """

    DEFAULT_ADDRESS = "http://www.ebi.ac.uk/arrayexpress/{format}/v2/"
    DEFAULT_FORMAT = "json"
    DEFAULT_CACHE = serverfiles.localpath(
        "ArrayExpress", "ArrayExpressCache.shelve")

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
        """Format the query arguments in `kwargs`.

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
            # TODO check if val contains a datetime.date object
            # assert proper format
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

        arg_items = sorted(
            arg_items,
            key=lambda arg: self._ARGS_ORDER.index(arg[0])
                            if arg[0] in self._ARGS_ORDER else 100
        )

        for key, value in arg_items:
            if key == "format":
                continue  # format is handled in query_url
            if key not in ARRAYEXPRESS_FIELDS:
                raise ValueError("Invalid argument name: '{0}'".format(key))
            if value is not None and value != []:
                fmt = formaters.get(key, format_default)
                value = fmt(value)
                parts.append("{0}={1}".format(key, value))

        return "&".join(parts)

    def query_url(self, what="experiments", **kwargs):
        """Return a formatted query URL for the query arguments.

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
        """Return query URL of formatted experiments for the query arguments.
        """
        return self.query_url("experiments", **kwargs)

    def query_url_files(self, **kwargs):
        """ Return query URL of formatted experiments for the query arguments.
        """
        return self.query_url("files", **kwargs)

    def query_experiment(self, **kwargs):
        """Return an open stream to the experiments query results. 
           Takes the same arguments as the :obj:`query_experiments` function.
        """
        url = self.query_url_experiments(**kwargs)
        stream = self._cache_urlopen(url, timeout=self.timeout)
        return stream

    def query_files(self, **kwargs):
        """Return an open stream to the files query results.
           Takes the same arguments as the :obj:`query_files` function.
        """
        url = self.query_url_files(**kwargs)
        stream = self._cache_urlopen(url, timeout=self.timeout)
        return stream

    def open_file(self, accession, kind="raw", ext=None):
        """ Return a file handle to experiment data.
        
        :param str accession:
        :param str kind: Experiment data type.
        
        Possible values for the parameter `kind`:
            - raw: return the raw data if available
            - processed: return the processed data if available
            - biosamples: a png or svg design image
            - idf: investigation description
            - adf: array design description
            - mageml: MAGE-ML file

        Example::

            >>> raw_file = conn.open_file("E-TABM-1087", kind="raw")
            >>> processed_file = conn.open_file("E-TABM-1087", kind="processed")

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
                return _fake_closing({})
        elif hasattr(self.cache, "close"):
            return closing(self.cache)
        elif self.cache is None:
            return _fake_closing({})
        else:
            return _fake_closing(self.cache)


@contextmanager
def _fake_closing(obj):
    yield obj


def query_experiments(keywords=None, accession=None, array=None, ef=None,
                      efv=None, expdesign=None, exptype=None,
                      gxa=None, pmid=None, sa=None, species=None,
                      expandefo=None, directsub=None, assaycount=None,
                      efcount=None, samplecount=None, rawcount=None,
                      fgemcount=None, miamescore=None, date=None,
                      format="json", wholewords=None, connection=None):
    """ Query Array Express experiments.

    :param keywords: A list of keywords to search (e.g. ``['gliobastoma']``).
    :param accession: Search by experiment accession (e.g. ``'E-MEXP-31'``).
    :param array: Search by array design name or accession (e.g. ``'A-AFFY-33'``).
    :param ef: Experimental factor (names of main variables of experiments).
    :param efv: Experimental factor value (Has EFO expansion).
    :param expdesign: Experiment design type (e.g. ``["dose", "response"]``).
    :param exptype: Experiment type (e.g. 'RNA-Seq', has EFO expansion).
    :param gxa: If True limit the results to the Gene Expression Atlas only.
    :param pmid: Search by PubMed identifier.
    :param sa: Sample attribute values (e.g. ``'fibroblast'``, has EFO expansion).
    :param species: Search by species (e.g. ``'Homo sapiens'``, has EFO expansion)
    :param expandefo: If True expand the search terms with all its child terms
        in the Experimental Factor Ontology (EFO_) (e.g. ``keywords="cancer"``
        will be expanded to include for synonyms and sub types of cancer).
    :param directsub: If True return only experiments submitted directly to
        Array Express; if False return only experiments imported from GEO
        database; if None (default) return both.
    :param assaycount: A two tuple (min, max) for filter on the number of
        assays (e.g. (1, 5) will return only experiments with at least 1
        and no more then 5 assays).
    :param efcount: Filter on the number of experimental factors (e.g. (1, 5)).
    :param sacount: Filter on the number of sample attribute categories.
    :param rawcount: Filter on the number or raw files.
    :param fgemcount: Filter on the number of final gene expression matrix
        (processed data) files.
    :param miamescore: Filter on the MIAME complience score (max 5).
    :param date: Filter by release date.

    >>> query_experiments(species="Homo sapiens", ef="organism_part", efv="liver")
    {u'experiments': ...

    .. _EFO: http://www.ebi.ac.uk/efo/

    """
    if connection is None:
        connection = ArrayExpressConnection()

    stream = connection.query_experiment(
        keywords=keywords, accession=accession,
        array=array, ef=ef, efv=efv, expdesign=expdesign, exptype=exptype,
        gxa=gxa, pmid=pmid, sa=sa, species=species, expandefo=expandefo,
        directsub=directsub, assaycount=assaycount, efcount=efcount,
        samplecount=samplecount, rawcount=rawcount, fgemcount=fgemcount,
        miamescore=miamescore, date=date, format=format,
        wholewords=wholewords
    )

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
    """Query Array Express files. See :obj:`query_experiments` for the
    arguments.

    >>> query_files(species="Mus musculus", ef="developmental_stage",
    ...             efv="embryo", format="xml")
    <xml.etree.ElementTree.ElementTree object ...

    """
    if connection is None:
        connection = ArrayExpressConnection()

    stream = connection.query_files(
        keywords=keywords, accession=accession,
        array=array, ef=ef, efv=efv, expdesign=expdesign, exptype=exptype,
        gxa=gxa, pmid=pmid, sa=sa, species=species, expandefo=expandefo,
        directsub=directsub, assaycount=assaycount, efcount=efcount,
        samplecount=samplecount, rawcount=rawcount, fgemcount=fgemcount,
        miamescore=miamescore, date=date, format=format,
        wholewords=wholewords
    )

    if format == "json":
        return parse_json(stream)
    else:
        return parse_xml(stream)


"""
MAGE-TAB convenience functions, classes
=======================================
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
    lines = [line.split("\t") for line in lines
             if line and not line.startswith("#")]
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
    lines = [line.split("\t") for line in lines
             if line.strip() and not line.startswith("#")]
    header = [h for h in lines[0] if h]
    rows = [line[:len(header)] for line in lines[1:]]
    assert(all([len(r) == len(header) for r in rows]))
    return header, rows


def parse_adf(file):
    raise NotImplementedError


def parse_data_matrix(file):
    """Parse the MAGE-TAB processed data matrix. Return a tuple where the
    elements are:
        - a (header REF, header values) tuple (e.g.
          ("Hybridization REF", ["Stage1", "Stage2", ...]) )
        - a list of quantitation type for header values (e.g.
          ["log2 ratio", "log2 ratio", ...])
        - a (row REF, row names list) tuple ("e.g.
          ("Reporter REF", ["Probe 1", "Probe 2", ...]) )
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
    r"""Investigation design (contains the contents of the .idf).

    >>> idf_file = StringIO(
    ...     'Investigation Title\tfoo investigation\n' +
    ...     'Experimental Design\tfubar\tsnafu\n' +
    ...     'SDRF File\tfoobar.sdrf\n'
    ... )
    >>> idf = InvestigationDesign(idf_file)
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
        try:
            return self[tag]
        except KeyError:
            pass

        ttag = self.transform_tag(tag)
        if hasattr(self, ttag):
            return getattr(self, ttag)
        else:
            raise KeyError(tag)


class SampleDataRelationship(object):
    """Sample-Data Relationship (contains the contents of the .sdrf file).

    """

    # Nodes an edges
    NODES_EDGES = [
        "Source Name", "Sample Name", "Extract Name",
        "Labeled Extract Name", "Hybridization Name",
        "Assay Name", "Scan Name", "Normalization Name",
        "Array Data File", "Derived Array Data File",
        "Array Data Matrix File", "Derived Array Data Matrix File",
        "Image File", "Protocol REF"]

    # Attributes for nodes and edges
    NODE_EDGE_ATTRIBUTES = {
        "Source Name": [
            "Characteristics", "Provider", "Material Type",
            "Description", "Comment"],
        "Sample Name": [
            "Characteristics", "Material Type", "Description", "Comment"],
        "Extract Name": [
            "Characteristics", "Material Type", "Description", "Comment"],
        "Labeled Extract Name": [
            "Characteristics", "Material Type", "Description",
            "Label", "Comment"],
        "Hybridization Name": [
            "Array Design File", "Array Design REF", "Comment"],
        "Assay Name": ["Technology Type", "Comment"],
        "Scan Name": ["Comment"],
        "Normalization Name": ["Comment"],
        "Array Data File": ["Comment"],
        "Derived Array Data File": ["Comment"],
        "Array Data Matrix File": ["Comment"],
        "Derived Array Data Matrix File": ["Comment"],
        "Image File": ["Comment"],
        "Protocol REF": [
            "Term Source REF", "Parameter", "Performer", "Date", "Comment"]
    }

    # Attributes
    ATTRIBUTE_COLUMNS = {
        "Characteristics []": ["Unit", "Term Source REF"],
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
        """ Transform the tag into a proper Python attribute name by
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
        """ Return the Source subsection.
        """
        return self._subsection("Source Name")

    def source_name(self):
        """ Return the Source Name subsection.
        """
        return self._column("Source Name")

    def sample(self):
        """ Return the Sample subsection.
        """
        return self._subsection("Sample Name")

    def sample_name(self):
        """ Return the Sample Name subsection
        """
        return self._column("Sample Name")

    def extract(self):
        """ Return the Extract subsection.
        """
        return self._subsection("Extract Name")

    def extract_name(self):
        """ Return the Extract Name subsection.
        """
        return self._column("Extract Name")

    def labeled_extract(self):
        """ Return the Labeled Extract subsection.
        """
        return self._subsection("Labeled Extract Name")

    def labeled_extract_name(self):
        """ Return the Labeled Extract Name subsection.
        """
        return self._column("Labeled Extract Name")

    def hybridization(self):
        """ Return the Hybridization subsection.
        """
        return self._subsection("Hybridization Name")

    def hybridization_name(self):
        """ Return the Hybridization Name subsection.
        """
        return self._column("Hybridization Name")

    def assay(self):
        """ Return the Assay subsection.
        """
        return self._subsection("Assay Name")

    def assay_name(self):
        """ Return the Assay Name subsection.
        """
        return self._column("Assay Name")

    def scan(self):
        """ Return the Scan subsection.
        """
        return self._subsection("Scan Name")

    def scan_name(self):
        """ Return the Scan name subsection.
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
        """ Return the Array Data subsection.
        """
        return self._subsection("Array Data File")

    def array_data_file(self):
        """ Return the Array Data File subsection.
        """
        return self._column("Array Data File")

    def derived_array_data(self):
        """ Return the Derived Array Data subsection.
        """
        return self._subsection("Derived Array Data File")

    def derived_array_data_file(self):
        """ Return the Derived Array Data File subsection.
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
    """ Array design (contains the contents of the .adf file).

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
    """ Are the strings in items continuous numbers. 
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
    """ Load a single processed matrix file into an :obj:`Orange.data.Table`.
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
    is_float = numpy.frompyfunc(_is_float, 1, 1)  # an numpy ufunc

    for header_name, quant, column in zip(header, quant_type, matrix.T):
        if _is_continuous(column):
            feature = Orange.feature.Continuous(header_name)
            # relace all non parsable floats with '?'
            column[numpy.where(1 - is_float(column))] = "?"
        else:
            values = set(column)
            feature = Orange.feature.Discrete(
                header_name, values=sorted(values)
            )
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
        pattern = re.compile(
            r"((Characteristics)|(Factor Value)|(Parameter Value)) " +
            r"\[(?P<name>.*?)\].*"
        )
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
    """Convert an MAGE-TAB annotated experiment into an 
    :obj:`Orange.data.Table` (assumes all the associated MAGE-TAB 
    files are in the same directory.

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
    """ Stack the tables horizontally.
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


class ArrayExpressExperiment(object):

    """
    An convenience class representing an Array Express Experiment.

    >>> ae = ArrayExpressExperiment("E-MEXP-2917")
    >>> print ae.name
    Characterization of Data Variation in Gene Expression Profiling of ...
    >>> for file in ae.files:
    ...     print file["name"], file["url"]
    E-MEXP-2917...
    >>> table = ae.fgem_to_table() # Retrieve the experiment data table

    :param str accession:
        The experiment accession id.

    """

    def __init__(self, accession, connection=None):
        #: Experiment accession id
        self.accession = accession
        #: A list of all species subjugated to torture in this experiment
        self.species = []
        #: True if Array Express provides raw data files.
        self.rawdatafiles = False
        #: True if Array Express provides processed data files.
        self.fgemdatafiles = False
        #: A list of all sample attributes.
        self.sampleattributes = []
        #: A list of experimental factors
        self.experimentalfactors = []

        self.connection = connection
        self._etree = tree = query_experiments(
            accession=accession, connection=self.connection, format="xml")
        experiments = tree.findall("experiment")
        # find the exact match (more then one experiment can be
        # listed in the query result)
        experiments = [e for e in experiments
                       if e.find("accession").text.strip() == accession]
        self._experiment = experiment = experiments[0]

        #: A list of all species subjugated to torture in this experiment
        self.species = [e.text for e in experiment.findall("species")]
        bool_values = {"true": True, "false": False}
        self.rawdatafiles = bool_values[experiment.find("rawdatafiles")
                                        .get("available", "false")]
        self.fgemdatafiles = bool_values[experiment.find("processeddatafiles")
                                         .get("available", "false")]

        for sa in experiment.findall("sampleattribute"):
            category = sa.find("category").text.strip()
            values = [val.text for val in sa.findall("value")]
            self.sampleattributes.append((category, values))

        for ef in experiment.findall("experimentalfactor"):
            name = ef.find("name").text.strip()
            values = [val.text.strip() for val in ef.findall("values")]
            self.experimentalfactors.append((name, values))

        self.miamescores = _dictify(experiment.find("miamescores"))

        self.id = experiment.find("id").text
        self.secondaryaccession = getattr(
            experiment.find("secondaryaccession"), "text", None
        )
        self.name = experiment.find("name").text
        self.experimenttype = experiment.find("experimenttype").text.strip()
        self.releasedate = experiment.find("releasedate").text
        self.lastupdatedate = getattr(
            experiment.find("lastupdatedate"), "text", None
        )
        self.samples = int(experiment.find("samples").text)
        self.assays = int(experiment.find("assays").text)

        self.arraydesign = \
            [_dictify(e) for e in experiment.findall("arraydesign")]

        self.bioassaydatagroups = \
            [_dictify(group)
             for group in experiment.findall("bioassaydatagroup")]
        self.bibliography = \
            [_dictify(e) for e in experiment.findall("bibliography")]
        self.provider = [_dictify(e) for e in experiment.findall("provider")]

        self.experimentdesign = []
        for expd in experiment.findall("experimentdesign"):
            self.experimentdesign.append(expd.text)

        self.description = \
            [_dictify(e) for e in experiment.findall("description")]

        tree = query_files(accession=self.accession,
                           format="xml",
                           connection=self.connection)
        experiments = tree.findall("experiment")
        experiments = [e for e in experiments
                       if e.find("accession").text.strip() == accession]
        experiment = experiments[0]
        files = experiment.findall("file")

        #: A list of file descriptions (dict instances) available provided
        #: by Array Express
        self.files = [_dictify(file) for file in files]

    def _download_file(self, url, extract=True):
        """ Download the `file` from the ArrayExpress into a local
        repository directory.

        """
        rest, basename = posixpath.split(url)
        dirname = posixpath.basename(rest)
        repo_dir = serverfiles.localpath("ArrayExpress", dirname)
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
        """ Is the `url` stored in the local repository?
        """
        return os.path.exists(self._local_filepath(url))

    def _local_filepath(self, url):
        """ Return the local file path for url.
        """
        rest, basename = posixpath.split(url)
        dirname = posixpath.basename(rest)
        return serverfiles.localpath(
                    "ArrayExpress", os.path.join(dirname, basename))

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
            extension_match = (extension == file.get("extension") or
                               extension is None)

            if kind_match and extension_match:
                res.append(file)
        return res

    def array_design(self):
        """ Return a list of :obj:`ArrayDesign` instances used in this experiment.
        """
        files = [f for f in self.files if f.get("kind") == "adf" and
                 f.get("extension") == "txt"]

        array_design = []
        for file in files:
            url = file.get("url")
            if not self._is_local(url):
                self._download_file(url)
            array_design.append(ArrayDesign(self._open(url)))
        return array_design

    def investigation_design(self):
        """ Return an :obj:`InvestigationDesign` for this experiment.
        """
        files = [f for f in self.files if f.get("kind") == "idf" and \
                 f.get("extension") == "txt"]
        if not files:
            raise ValueError("The experiment '{0}' does not have an "
                             "investigation design file"
                             .format(self.accession))
        file = files[0]
        return InvestigationDesign(self._open(file.get("url")))

    def sample_data_relationship(self):
        """
        Return an :obj:`SampleDataRelationship` instance describing this experiment.
        """
        files = [f for f in self.files if f.get("kind") == "sdrf" and \
                 f.get("extension") == "txt"]
        if not files:
            raise ValueError("The experiment '{0}' does not have an sample "
                             "and data relationship file"
                             .format(self.accession))
        file = files[0]
        return SampleDataRelationship(self._open(file.get("url")))

    def fgem_to_table(self):
        """ Retrieve the processed matrix from the Array Express FTP
        server and convert it to a :class:`Orange.data.Table`.

        """
        assert(self.fgemdatafiles)
        repo_dir = serverfiles.localpath("ArrayExpress", self.accession)
        # Find the file listing the data matrix files
        # (should be in sdrf but sometimes it is in 2column file only, why?)
        sdrf = self._search_files("sdrf", "txt")
        if sdrf:
            sdrf = SampleDataRelationship(self._open(sdrf[0].get("url")))
            if "Derived Array Data Matrix File" not in sdrf.header:
                twocol = self._search_files("twocolumn", "txt")
                if twocol:
                    sdrf = SampleDataRelationship(
                        self._open(twocol[0].get("url"))
                    )
        matrix_file = self._search_files("fgem")[0]
        self._open(matrix_file.get("url"))

        idf_file = self._search_files("idf", "txt")[0]
        self._open(idf_file.get("url"))  # To download if not cached
        return mage_tab_to_orange(os.path.join(repo_dir, idf_file.get("name")))


def test():
    conn = ArrayExpressConnection()
    import doctest
    foo = type("foo", (object,), {})()
    bar = type("bar", (object,), {})()

    doctest.testmod(optionflags=doctest.ELLIPSIS,
                    extraglobs={"conn": conn, "foo": foo, "bar": bar})


if __name__ == "__main__":
    test()
