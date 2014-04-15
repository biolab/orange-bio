import os
import urllib2
import anydbm
import shelve
import itertools
import warnings

from functools import wraps
from collections import namedtuple
from operator import itemgetter
from contextlib import closing
from xml.dom import pulldom
from StringIO import StringIO

from Orange.utils import environ


class BioMartError(Exception):
    pass


class BioMartServerError(BioMartError):
    pass


class DatabaseNameConflictError(BioMartError):
    pass


class BioMartQueryError(BioMartError):
    pass


_BaseAttribute = namedtuple(
    "Attribute",
    ["internal_name", "name", "description", "4", "format", "table", "7"],
    rename=True
)


class Attribute(_BaseAttribute):
    """
    An attribute in a BioMart data-set.
    """
    @classmethod
    def from_string(cls, string):
        return cls(string.split("\t"))

    # Re-define the getters for documentation
    internal_name = property(
        itemgetter(0), doc="Attribute's internal name")

    name = property(
        itemgetter(1), doc="Human readable name.")

    description = property(
        itemgetter(2), doc="Human readable description")

    format = property(
        itemgetter(4), doc="Attribute format")

    # For back compatibility.
    internalName = internal_name


_BaseFilter = namedtuple(
    "Filter",
    ["internalName", "name", "values", "description", "type",
     "qualifiers", "6", "attribute_table"],
    rename=True
)


class Filter(_BaseFilter):
    """
    A filter on a BioMart data-set.
    """
    @classmethod
    def from_string(cls, string):
        return cls(string.split("\t"))

    internal_name = property(
        itemgetter(0), doc="Internal name")

    name = property(
        itemgetter(1), doc="Filter name")

    values = property(
        itemgetter(2), doc="List of possible filter values")

    description = property(
        itemgetter(3), doc="Filter description")


class XMLNode(object):

    def __init__(self, tag, attributes, children=[], data=""):
        self.tag = tag
        self.attributes = attributes
        self.children = children if children else []
        self.data = data if data else ""

    @classmethod
    def fromDOMNode(cls, node):
        return XMLNode(node.tagName, dict(node.attributes.items()))

    def _match(self, tag=None, **kwargs):
        match = True
        if tag is not None and self.tag != tag:
            match = False

        return match and all(self.attributes.get(key, None) == value
                             for key, value in kwargs.iteritems())

    def elements(self, tag=None, **kwargs):
        if self._match(tag, **kwargs):
            yield self

        for child in self.children:
            for el in child.elements(tag, **kwargs):
                yield el

    def elements_top(self, tag=None, **kwargs):
        """ Return only the top elements (do not return matching children
        of an already matching element)."
        """
        if self._match(tag, **kwargs):
            yield self
        else:
            for child in self.children:
                for el in child.elements_top(tag, **kwargs):
                    yield el

    def subelements(self, *args, **kwargs):
        return itertools.chain(*[c.elements(*args, **kwargs)
                                 for c in self.children])

    def subelements_top(self, *args, **kwargs):
        return itertools.chain(*[c.elements_top(*args, **kwargs)
                                 for c in self.children])

    def __iter__(self):
        """ Iterate over all sub elements.
        """
        return itertools.chain(*[iter(c) for c in self.children])

    def __repr__(self):
        return 'XMLNode("%s", %s)' % (self.tag, repr(self.attributes))


def parseXML(stream, parser=None):
    if isinstance(stream, basestring):
        events = pulldom.parseString(stream, parser)
    else:
        events = pulldom.parse(stream, parser)

    document = None
    chain = []
    for event, node in events:
        if event == "START_DOCUMENT":
            chain.append(XMLNode("DOCUMENT", {}))

        elif event == "START_ELEMENT":
            node = XMLNode.fromDOMNode(node)
            if chain:
                chain[-1].children.append(node)
            chain.append(node)

        elif event == "END_ELEMENT":
            chain.pop(-1)

        elif event == "CHARACTERS":
            chain[-1].data += node.data

        elif event == "END_DOCUMENT":
            document = chain.pop(-1)
    return document or chain[0]


def de_tab(text, sep="\t"):
    return [line.split(sep) for line in text.splitlines() if line.strip()]


def cached(func):
    @wraps(func)
    def f(self):
        if getattr(self, "_cache_" + func.__name__, None) is None:
            setattr(self, "_cache_" + func.__name__, func(self))
        return getattr(self, "_cache_" + func.__name__)
    return f


def drop_iter(generator):
    """Drop all elements from the iterator that raise an exception.
    """
    def _iter(generator):
        while True:
            try:
                yield generator.next()
            except StopIteration:
                raise StopIteration
            except Exception, ex:
                warnings.warn("An error occured during iteration:\n%s" %
                              str(ex), UserWarning, stacklevel=2)
    return list(_iter(generator))


DEFAULT_ADDRESS = "http://www.biomart.org/biomart/martservice"

DATA_CACHE = os.path.join(environ.buffer_dir, "biomart-data.cache.db")
META_CACHE = os.path.join(environ.buffer_dir, "biomart-mata.cache.db")


def checkBioMartServerError(response):
    if response.strip().startswith("Mart name conflict"):
        raise DatabaseNameConflictError(response)
    elif response.strip().startswith("Problem retrieving datasets"):
        raise BioMartError(response)
    elif response.startswith("non-BioMart die():"):
        raise BioMartServerError(response)


class BioMartConnection(object):

    """ A connection to a BioMart martservice server.

    >>> connection = BioMartConnection(
    ...     "http://www.biomart.org/biomart/martservice")
    >>> response = connection.registry()
    >>> response = connection.datasets(mart="ensembl")

    """
    FOLLOW_REDIRECTS = False

    def __init__(self, address=None):

        self.address = address if address is not None else DEFAULT_ADDRESS

        self._error_cache = {}

    def _get_cache_for_args(self, kwargs):
        if "type" in kwargs:
            return self._open_meta_cache
        elif "query" in kwargs:
            return self._open_data_cache
        else:
            raise ValueError

    def _open_data_cache(self, flag="r"):
        return self._open_cache(DATA_CACHE, flag=flag)

    def _open_meta_cache(self, flag="r"):
        return self._open_cache(META_CACHE, flag=flag)

    def _open_cache(self, filename, flag="r"):
        # Ensure the cache db actualy exists (r,w flags do not create it)
        if flag in ["r", "w"]:
            try:
                return shelve.open(filename, flag)
            except anydbm.error:
                shelve.open(filename, "c").close()
                return shelve.open(filename, flag)
        else:
            return shelve.open(filename, flag)

    def request_url(self, **kwargs):
        order = ["type", "dataset", "mart", "virtualSchema", "query"]
        items = sorted(
           kwargs.iteritems(),
           key=lambda item: order.index(item[0]) if item[0] in order else 10
        )
        query = "&".join("%s=%s" % (key, urllib2.quote(value))
                         for key, value in items if key != "POST")

        return self.address + "?" + query
#         return url.replace(" ", "%20")

    def request(self, **kwargs):
        url = self.request_url(**kwargs)
        cache_key = url
        if isinstance(url, unicode):
            cache_key = url.encode("utf-8")

        if cache_key in self._error_cache:
            raise self._error_cache[cache_key]

        cache_open = self._get_cache_for_args(kwargs)

        with closing(cache_open(flag="r")) as cache:
            response_cached = cache_key in cache
            response = cache[cache_key] if response_cached else None

        if not response_cached:
            try:
                response = urllib2.urlopen(url).read()
                checkBioMartServerError(response)
            except Exception, ex:
                self._error_cache[url] = ex
                raise ex

            with closing(cache_open(flag="w")) as cache:
                cache[cache_key] = response

        return StringIO(response)

    def registry(self, **kwargs):
        return self.request(type="registry")

    def datasets(self, mart="ensembl", **kwargs):
        return self.request(type="datasets", mart=mart, **kwargs)

    def attributes(self, dataset="oanatinus_gene_ensembl", **kwargs):
        return self.request(type="attributes", dataset=dataset, **kwargs)

    def filters(self, dataset="oanatinus_gene_ensembl", **kwargs):
        return self.request(type="filters", dataset=dataset, **kwargs)

    def configuration(self, dataset="oanatinus_gene_ensembl", **kwargs):
        return self.request(type="configuration", dataset=dataset, **kwargs)

    def clear_cache(self):
        with closing(self._open_data_cache(flag="n")):
            pass

        with closing(self._open_meta_cache(flag="n")):
            pass

        self._error_cache.clear()

    # Back compatibility
    clearCache = clear_cache


class BioMartRegistry(object):

    """ A class representing a BioMart registry. Arguments:

    :param stream: A file like object with xml registry or a
          BioMartConnection instance

    >>> registry = tioMartRegistry(connection)
    >>> for schema in registry.virtual_schemas():
    ...    print schema.name
    ...
    default

    """

    def __init__(self, stream):
        self.connection = stream if isinstance(stream, BioMartConnection) \
                          else None
        if self.connection:
            if hasattr(self.connection, "_registry"):
                self.__dict__ = self.connection._registry.__dict__
                return
            else:
                stream = self.connection.registry()
                self.connection._registry = self
        else:
            stream = open(stream, "rb") if isinstance(stream, basestring) \
                     else stream
        self.registry = self.parse(stream)

    @cached
    def virtual_schemas(self):
        """Return a list of :obj:`BioMartVirtualSchema` instances representing
        each schema.
        """
        schemas = [schema.attributes.get(name, "default")
                   for name in self.registry.elements("virtualSchema")]
        if not schemas:
            schemas = [BioMartVirtualSchema(self.registry, name="default",
                                            connection=self.connection)]
        else:
            schemas = [BioMartVirtualSchema(
                            schema,
                            name=schema.attributes.get("name", "default"),
                            connection=self.connection)
                       for schema in self.registry.elements("virtualSchema")]
        return schemas

    def virtual_schema(self, name):
        """ Return a named virtual schema.
        """
        for schema in self.virtual_schemas():
            if schema.name == name:
                return schema
        raise ValueError("Unknown schema name '%s'" % name)

    @cached
    def marts(self):
        """Return a list off all 'mart' instances (:obj:`BioMartDatabase`)
        regardless of their virtual schemas.
        """
        return reduce(list.__add__,
                      drop_iter(schema.marts()
                                for schema in self.virtual_schemas()),
                      [])

    def mart(self, name):
        """ Return a named mart.
        """
        for mart in self.marts():
            if mart.name == name:
                return mart
        raise ValueError("Unknown mart name '%s'" % name)

    def databases(self):
        """Same as :obj:`marts`.
        """
        return self.marts()

    def datasets(self):
        """Return a list of all datasets (:obj:`BioMartDataset`) from all
        marts regardless of their virtual schemas.
        """
        return reduce(list.__add__,
                      drop_iter(mart.datasets()
                                for mart in self.marts()),
                      [])

    def dataset(self, internalName, virtualSchema=None):
        """ Return a BioMartDataset instance that matches the `internalName`.
        """
        if virtualSchema is not None:
            schema = self.virtual_schema(virtualSchema)
            return schema.dataset(internalName)
        else:
            for d in self.datasets():
                if d.internalName == internalName:
                    return d
            raise ValueError("Unknown dataset name: '%s'" % internalName)

    def query(self, **kwargs):
        """Return an initialized :obj:`BioMartQuery` with registry set to self.
        Pass additional arguments to BioMartQuery.__init__ with keyword
        arguments.
        """
        return BioMartQuery(self, **kwargs)

    def links_between(self, exporting, importing, virtualSchema="default"):
        """ Return all links between `exporting` and `importing` datasets in the
        `virtualSchema`.
        """
        schema = [schema for schema in self.virtual_schemas()
                  if schema.name == virtualSchema][0]
        return schema.links_between(exporting, importing)

    def get_path(self, exporting, importing):
        raise NotImplementedError

    def __iter__(self):
        return iter(self.marts())

    def _find_mart(self, name):
        try:
            return [mart for mart in self.marts()
                    if name in [getattr(mart, "name", None),
                                getattr(mart, "internalName", None)]][0]
        except IndexError, ex:
            raise ValueError(name)

    def __getitem__(self, name):
        return self._find_mart(name)

    def search(self, string, relevance=False):
        results = []
        for mart in self.marts():
            results.extend([(rel, mart, dataset)
                            for rel, dataset in mart.search(string, True)])
        results = sorted(results, reverse=True)
        if not relevance:
            results = [(mart, dataset) for _, mart, dataset in results]
        return results

    @classmethod
    def parse(cls, stream, parser=None):
        """Parse the registry file like object and return a DOM like
        description (:obj:`XMLNode`).
        """
        xml = stream.read()
        doc = parseXML(xml, parser)
        return doc

parseRegistry = BioMartRegistry.parse


class BioMartVirtualSchema(object):
    """ A virtual schema.
    """

    def __init__(self, locations=None, name="default", connection=None):
        self.locations = locations
        self.name = name
        self.connection = connection

    @cached
    def marts(self):
        """Return a list off all 'mart' instances (:obj:`BioMartDatabase`)
        in this schema.
        """
        return drop_iter(
            BioMartDatabase(connection=self.connection,
                            **dict((str(key), val)
                                   for key, val in loc.attributes.items()))
            for loc in self.locations.elements("MartURLLocation")
        )

    def databases(self):
        """Same as :obj:`marts`.
        """
        return self.marts()

    @cached
    def datasets(self):
        """Return a list of all datasets (:obj:`BioMartDataset`)
        from all marts in this schema.
        """
        return reduce(list.__add__,
                      drop_iter(mart.datasets()  for mart in self.marts()),
                      [])

    def dataset(self, internalName):
        """Return a dataset with matching `internalName`.
        """
        for dataset in self.datasets():
            if dataset.internalName == internalName:
                return dataset
        raise ValueError("Unknown data set name '%s'" % internalName)

    def links_between(self, exporting, importing):
        """Return a list of link names from `exporting` dataset to
        `importing` dataset.
        """
        exporting = self[exporting]
        importing = self[importing]
        exporting = exporting.configuration().exportables()
        importing = importing.configuration().importables()
        exporting = set([(ex.linkName, getattr(ex, "linkVersion", ""))
                         for ex in exporting])
        importing = set([(im.linkName, getattr(ex, "linkVersion", ""))
                         for im in importing])
        links = exporting.intersection(importing)
        return [link[0] for link in links]

    def links(self):
        """Return a list of (linkName, linkVersion) tuples defined by
        datasets in the schema.
        """
        pass

    def query(self, **kwargs):
        """ Return an initialized :obj:`BioMartQuery` with registry
        and `virtualSchema` set to self. Pass additional arguments to
        :obj:`BioMartQuery.__init__` with keyword arguments

        """
        return BioMartQuery(self, virtualSchema=self.name, **kwargs)

    def __iter__(self):
        return iter(self.marts())

    def __getitem__(self, key):
        try:
            return self.dataset(key)
        except ValueError:
            raise KeyError(key)


class BioMartDatabase(object):
    """A BioMart 'mart' instance.

    :param str name: Name of the mart instance.
    :param str virtualSchema: Name of the virtualSchema of the dataset.
    :param BioMartConnection: An optional BioMartConnection instance.

    """
    host = "www.biomart.org"
    path = "/biomart/martservice"

    def __init__(self, name="ensembl", virtualSchema="default",
                 connection=None, database="ensembl_mart_60", default="1",
                 displayName="ENSEMBL GENES 60 (SANGER UK)",
                 host="www.biomart.org", includeDatasets="", martUser="",
                 path="/biomart/martservice", port="80",
                 serverVirtualSchema="default", visible="1", **kwargs):

        self.name = name
        self.virtualSchema = virtualSchema
        self.database = database,
        self.default = default
        self.displayName = displayName
        self.host = host
        self.includeDatasets = includeDatasets
        self.martUser = martUser
        self.path = path
        self.port = port
        self.serverVirtualSchema = serverVirtualSchema
        self.visible = visible
        self.__dict__.update(kwargs.items())  # why?

        if connection is None:
            connection = BioMartConnection()

        if kwargs.get("redirect", None) == "1" and \
                BioMartConnection.FOLLOW_REDIRECTS:
            redirect = BioMartConnection(
                "http://" + self.host + ":" + self.port + self.path,
                cache=connection.cache)
            try:
                registry = redirect.registry()
                connection = redirect
            except urllib2.HTTPError, ex:
                warnings.warn("'%s' is not responding!, using the default "
                              "original connection. %s" %
                              (redirect.address, str(ex)))

        self.connection = connection

    @cached
    def _datasets_index(self):
        keys = ["datasetType", "internalName", "displayName", "visible",
                "assembly", "_5", "_6", "virtualSchema", "date"]

        try:
            datasets = self.connection.datasets(
                mart=self.name, virtualSchema=self.virtualSchema).read()
        except BioMartError, ex:
            if self.virtualSchema == "default":
                datasets = self.connection.datasets(mart=self.name).read()
            else:
                raise
        datasets = de_tab(datasets)
        return [dict(zip(keys, line)) for line in datasets]

    @cached
    def datasets(self):
        """Return a list of all datasets (:obj:`BioMartDataset`)
        in this database.
        """
        return drop_iter(BioMartDataset(mart=self.name,
                                        connection=self.connection, **dataset)
                         for dataset in self._datasets_index())

    def dataset_attributes(self, dataset, **kwargs):
        """ Return a list of dataset attributes.
        """
        dataset = self._find_dataset(dataset)
        return dataset.attributes()

    def dataset_filters(self, dataset, **kwargs):
        """ Return a list of dataset filters.
        """
        dataset = self._find_dataset(dataset)
        return dataset.filters()

    def dataset_query(self, dataset, filters=[], attributes=[]):
        """ Return an dataset query based on dataset, filters and attributes.
        """
        dataset = self._find_dataset(dataset)
        return BioMartQuery(self.con, [dataset], filters, attributes).run()

    def _find_dataset(self, internalName):
        for dataset in self.datasets():
            if dataset.internalName == internalName:
                return dataset
        raise ValueError("Uknown dataset name '%s'" % internalName)

    def __getitem__(self, name):
        try:
            return self._find_dataset(name)
        except ValueError:
            raise KeyError(name)

    def search(self, string, relevance=False):
        strings = string.lower().split()
        results = []
        for dataset, conf in drop_iter((dataset, dataset.configuration())
                                       for dataset in self.datasets()):
            trees = (conf.attribute_pages() +
                     conf.attribute_groups() +
                     conf.attribute_collections() +
                     conf.attributes())

            names = (" ".join([getattr(tree, "displayName", "")
                               for tree in trees])
                     .lower())
            count = sum([names.count(s) for s in strings])
            if count:
                results.append((count, dataset))
        return results


class BioMartDataset(object):
    """
    A BioMart dataset (returned by :obj:`BioMartDatabase`).
    """

    # these attributes store the result of
    # ".../martservice?type=datasets&mart=..." query
    FIELDS = ["datasetType", "internalName", "displayName", "visible",
              "assembly", "_5", "_6", "virtualSchema", "date"]

    def __init__(self, mart="ensembl", internalName="hsapiens_gene_ensembl",
                 virtualSchema="default", connection=None,
                 datasetType="TableSet", displayName="", visible="1",
                 assembly="", date="", **kwargs):
        self.connection = \
            connection if connection is not None else BioMartConnection()
        self.mart = mart
        self.internalName = internalName
        self.virtualSchema = virtualSchema
        self.datasetType = datasetType
        self.displayName = displayName
        self.visible = visible
        self.assembly = assembly
        self.date = date
        self.__dict__.update(kwargs)
        self._attributes = None
        self._filters = None

    @cached
    def attributes(self):
        """Return a list of available attributes for this dataset
        (:obj:`Attribute`).
        """
        stream = self.connection.attributes(
            dataset=self.internalName, virtualSchema=self.virtualSchema
        )
        response = stream.read()
        lines = response.splitlines()
        return [Attribute.from_string(line) for line in lines if line.strip()]

    @cached
    def filters(self):
        """Return a list of available filters for this dataset
        (:obj:`Filter`).
        """
        stream = self.connection.filters(
            dataset=self.internalName, virtualSchema=self.virtualSchema
        )
        response = stream.read()
        lines = response.splitlines()
        return [Filter.from_string(line) for line in lines if line.strip()]

    @cached
    def configuration(self, parser=None):
        """Return the configuration tree for this dataset
        (:obj:`DatasetConfig`).
        """
        stream = self.connection.configuration(
            dataset=self.internalName, virtualSchema=self.virtualSchema
        )
        response = stream.read()
        doc = parseXML(response, parser)
        config = list(doc.elements("DatasetConfig"))[0]
        return DatasetConfig(BioMartRegistry(self.connection), config.tag,
                             config.attributes, config.children)

    def get_data(self, attributes=[], filters=[], unique=False):
        """Construct and run a :obj:`BioMartQuery` and return its results.
        """
        query = BioMartQuery(self.connection, dataset=self,
                             attributes=attributes, filters=filters,
                             uniqueRows=unique,
                             virtualSchema=self.virtualSchema)
        return query.run()

    def count(self, filters=[], unique=False):
        """Construct and run a :obj:`BioMartQuery` and count the number of returned
        lines.
        """
        query = BioMartQuery(self.connection, dataset=self, filters=filters,
                             uniqueRows=unique,
                             virtualSchema=self.virtualSchema)
        return query.get_count()

    def get_example_table(self, attributes=[], filters=[], unique=False):
        query = BioMartQuery(self.connection, dataset=self,
                             attributes=attributes, filters=filters,
                             uniqueRows=unique,
                             virtualSchema=self.virtualSchema)
        return query.get_example_table()


class BioMartQuery(object):
    """ Construct a query to run on a BioMart server.

    >>> BioMartQuery(connection,
    ...              dataset="hsapiens_gene_ensembl",
    ...              attributes=["ensembl_transcript_id",
    ...                          "chromosome_name"],
    ...              filters=[("chromosome_name", ["22"])]).get_count()
    1221
    >>> # Equivalent to
    >>> query = BioMartQuery(connection)
    >>> query.set_dataset("hsapiens_gene_ensembl")
    >>> query.add_filter("chromosome_name", "22")
    >>> query.add_attribute("ensembl_transcript_id")
    >>> query.add_attribute("chromosome_name")
    >>> query.get_count()
    1221

    """
    class XMLQuery(object):
        XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName = "%(virtualSchemaName)s" formatter = "%(formatter)s" header = "%(header)s" uniqueRows = "%(uniqueRows)s" count = "%(count)s" datasetConfigVersion = "%(version)s" >  
%(datasets_xml)s
</Query>"""
        version = "0.6"

        def __init__(self, query):
            self.query = query
            self._query = query._query
            self.__dict__.update(query.__dict__)

        def get_xml(self, count=None, header=False):
            datasets_xml = "\n".join(self.query_dataset(*query)
                                     for query in self._query)

            links_xml = self.query_links(self._query)
            datasets_xml += "\n" + links_xml

            count = self.count if count is None else count
            args = dict(datasets_xml=datasets_xml,
                        uniqueRows="1" if self.uniqueRows else "0",
                        count="1" if count else "",
                        header="1" if header else "0",
                        virtualSchemaName=self.virtualSchema,
                        formatter=self.format,
                        version=self.version)
            xml = self.XML % args
            return xml

        def query_attributes(self, attributes):
            def name(attr):
                if isinstance(attr, Attribute):
                    return attr.internalName
                elif isinstance(attr, AttributeDescription):
                    return attr.internalName
                else:
                    return str(attr)

            xml = "\n".join('\t\t<Attribute name = "%s" />' %
                            name(attr).lower()
                            for attr in attributes)
            return xml

        def query_filter(self, filters):
            def name(filter):
                if isinstance(filter, Filter):
                    return filter.internalName
                elif isinstance(filter, FilterDescription):
                    return getattr(filter, "field",
                                   getattr(filter, "internalName"))
                else:
                    return str(filter)

            def args(args):
                if isinstance(args, list):
                    return 'value="%s"' % ",".join(args)
                elif isinstance(args, basestring):
                    return 'value="%s"' % args
                elif isinstance(args, dict):
                    return " ".join(['%s="%s"' % (key, value)
                                     for key, value in args.iteritems()])

            xml = "\n".join(
                '\t\t<Filter name = "%s" %s />' % (name(fil), args(val))
                for fil, val in filters
            )
            return xml

        def query_dataset(self, dataset, attributes, filters):
            name, interface = (dataset.internalName, "default") if \
                isinstance(dataset, BioMartDataset) else (dataset, "default")
            xml = ('\t<Dataset name = "%s" interface = "%s" >\n' %
                   (name, interface))
            xml += self.query_attributes(attributes) + "\n"
            xml += self.query_filter(filters)
            xml += '\n\t</Dataset>'
            return xml

        def query_links(self, query):
            def name(dataset):
                if isinstance(dataset, BioMartDataset):
                    return getattr(dataset, "internalName",
                                   getattr(dataset, "name", None))
                else:
                    return str(dataset)

            if len(query) == 2:
                return ('\t<Links source="%s" target="%s"/>' %
                        (name(query[0][0]), name(query[1][0])))
            else:
                return ""

    class XMLQueryOld(XMLQuery):
        version = "0.4"

        def query_filter(self, filters):
            def name(filter):
                if isinstance(filter, Filter):
                    return filter.internalName
                elif isinstance(filter, FilterDescription):
                    return getattr(filter, "field",
                                   getattr(filter, "internalName"))
                else:
                    return str(filter)

            def value(value):
                return (str(value) if not isinstance(value, list)
                        else ",".join([str(v) for v in value]))

            xml = "\n".join(
                ('\t\t<ValueFilter name = "%s" value="%s" />' %
                 (name(fil), value(val)))
                for fil, val in filters)
            return xml

    def __init__(self, registry, virtualSchema="default", dataset=None,
                 attributes=[], filters=[], count=False, uniqueRows=False,
                 format="TSV"):
        if isinstance(registry, BioMartConnection):
            self.registry = BioMartRegistry(registry)
            self.virtualSchema = virtualSchema
        elif isinstance(registry, BioMartVirtualSchema):
            self.registry = registry.connection
            self.virtualSchema = registry.name
        else:
            self.registry = registry
            self.virtualSchema = virtualSchema

        self._query = []
        if dataset:
            self.set_dataset(dataset)
            for attr in attributes:
                self.add_attribute(attr)
            for filter, value in filters:
                self.add_filter(filter, value)
        self.count = count
        self.uniqueRows = uniqueRows
        self.format = format

    def set_dataset(self, dataset):
        self._query.append((dataset, [], []))

    def add_filter(self, filter, value):
        self._query[-1][2].append((filter, value))

    def add_attribute(self, attribute):
        self._query[-1][1].append(attribute)

    def get_count(self):
        count = self.run(count=True)
        if count.strip():
            count = int(count.strip())
        else:
            count = 0
        return count

    def run(self, count=None, header=False):
        query = (self.xml_query(count=count, header=header)
                     .replace("\n", "").replace("\t", ""))
        stream = self.registry.connection.request(query=query)
        stream = stream.read()
        if stream.startswith("Query ERROR:"):
            raise BioMartQueryError(stream)
        return stream

    def set_unique(self, unique=False):
        self.uniqueRows = unique

    def xml_query(self, count=None, header=False):
        self.version = version = "0.4"

        schema = self.virtualSchema

        dataset = self._query[0][0]
        dataset = dataset.internalName if isinstance(dataset, BioMartDataset) \
                  else dataset

        dataset = self.registry.dataset(
            dataset, virtualSchema=self.virtualSchema)

        conf = dataset.configuration()
        self.version = version = getattr(conf, "softwareVersion", "0.4")

        if version > "0.4":
            xml = self.XMLQuery(self).get_xml(count, header)
        else:
            xml = self.XMLQueryOld(self).get_xml(count, header)
        return xml

    def get_example_table(self):
        import orange
        data = self.run(count=False, header=True)

        if self.format.lower() == "tsv":
            header, data = data.split("\n", 1)
            domain = orange.Domain([orange.StringVariable(name)
                                    for name in header.split("\t")], False)
            data = [line.split("\t")
                    for line in data.split("\n") if line.strip()]
            return orange.ExampleTable(domain, data) if data else None
        elif self.format.lower() == "fasta":
            domain = orange.Domain(
                [orange.StringVariable("id"),
                 orange.StringVariable("sequence")],
                False)  # TODO: meaningful id
            examples = []
            from StringIO import StringIO
            from Bio import SeqIO
            for seq in SeqIO.parse(StringIO(data), "fasta"):
                examples.append([seq.id, str(seq.seq)])
            return orange.ExampleTable(domain, examples)
        else:
            raise BioMartError("Unsupported format: %" % self.format)


def get_pointed(self):
    if self.is_pointer():
        pointerDataset = self.pointerDataset
        if hasattr(self, "pointerAttribute"):
            name, getter = self.pointerAttribute, "AttributeDescription"
        elif hasattr(self, "pointerFilter"):
            name, getter = self.pointerFilter, "FilterDescription"

        dataset = self.registry.dataset(pointerDataset)

        conf = dataset.configuration()
        desc = list(conf.elements(getter, internalName=name))
        if desc:
            return dataset, desc[0]
        else:
            warnings.warn("Could not resolve pointer '%s' in '%s'" %
                          (name, pointerDataset), UserWarning, stacklevel=2)
            return None, None


def is_pointer(self):
    return hasattr(self, "pointerAttribute") or hasattr(self, "pointerFilter")


class ConfigurationNode(XMLNode):

    def __init__(self, registry, *args, **kwargs):
        XMLNode.__init__(self, *args, **kwargs)
        self.registry = registry
        self.__dict__.update([(self._name_mangle(name), value) \
                              for name, value in self.attributes.iteritems()])

        self.children = [self._factory(child) for child in self.children]

    def _name_mangle(self, name):
        if name in self.__dict__:
            return name + "_"
        else:
            return name

    def _factory(self, node):
        tag = node.tag
        class_ = globals().get(tag, None)
        if not class_:
            raise ValueError("Unknown node '%s;" % tag)
        else:
            return class_(self.registry, tag, node.attributes, node.children)

    def elements(self, tag=None, **kwargs):
        if isinstance(tag, type):
            tag = tag.__name__
        return XMLNode.elements(self, tag, **kwargs)

    def elements_top(self, tag=None, **kwargs):
        if isinstance(tag, type):
            tag = tag.__name__
        return XMLNode.elements_top(self, tag, **kwargs)

    def __repr__(self):
        return ('%s("%s", %s)' %
                (self.__class__.__name__, self.tag, repr(self.attributes)))


class DatasetConfig(ConfigurationNode):
    pass


class Exportable(ConfigurationNode):
    pass


class Importable(ConfigurationNode):
    pass


class FilterPage(ConfigurationNode):
    pass


class FilterGroup(ConfigurationNode):
    pass


class FilterCollection(ConfigurationNode):
    pass


class FilterDescription(ConfigurationNode):

    is_pointer = is_pointer
    get_pointed = get_pointed


class AttributePage(ConfigurationNode):
    pass


class AttributeGroup(ConfigurationNode):
    pass


class AttributeCollection(ConfigurationNode):
    pass


class AttributeDescription(ConfigurationNode):

    is_pointer = is_pointer
    get_pointed = get_pointed


class Option(ConfigurationNode):
    pass


class PushAction(ConfigurationNode):
    pass


class MainTable(ConfigurationNode):
    pass


class Key(ConfigurationNode):
    pass


if __name__ == "__main__":
    con = BioMartConnection("http://www.biomart.org/biomart/martservice")
    registry = BioMartRegistry(con)
    for schema in registry.virtual_schemas():
        print "Virtual schema '%s'" % schema.name
        for mart in schema.databases():
            print "\tMart: '%s' ('%s'):" % (mart.name, mart.displayName)
            for dataset in mart.datasets():
                print "\t\t Dataset '%s' %s' '%s'" % \
                      (dataset.datasetType, dataset.internalName,
                       dataset.displayName)

    database = BioMartDatabase(name="dicty", connection=con)
    datasets = database.datasets()
    print datasets
    dataset = datasets[2]
    configuration = dataset.configuration()
    attr = dataset.attributes()
    print attr
    filters = dataset.filters()
    print filters
    reg = BioMartRegistry(con)

    dataset = reg.dataset("scerevisiae_gene_ensembl")
    query = BioMartQuery(
         con, dataset="scerevisiae_gene_ensembl",
         attributes=["ensembl_transcript_id", "transcript_exon_intron"],
         filters=[("chromosome_name", "I"),
                  ("with_wikigene", {"excluded": "1"})],
         format="FASTA")
    print query.xml_query()
    print query.run()

    data = query.get_example_table()
    data.save("seq.tab")

    import doctest
    doctest.testmod(extraglobs={"connection": con, "registry": registry},
                    optionflags=doctest.ELLIPSIS)
