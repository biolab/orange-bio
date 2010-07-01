""" A python module for accessing BioMart services
Example::
    >>> connection = BioMartConnection("http://www.biomart.org/biomart/martservice")
    >>> reg = BioMartRegistry(connection)
    >>> for mart in reg.marts():
    ...    print mart.name, mart.datasets()
    ...
    ensembl [<obiBioMart.BioMartDataset ...
    >>> dataset = BioMartDataset(mart="ensembl", name="hsapiens_gene_ensembl", virtualSchema="default")
    >>> for attr in dataset.attributes():
    ...    print attr.name
    ...
    Ensembl Gene ID
    ...
    Constitutive Exon
    >>> data = dataset.get_data(attributes=["ensembl_gene_id", "ensembl_peptide_id"], filters=[("chromosome_name", "1")])
    
    >>> query = BioMartQuery(reg.connection, virtualSchema="default")
    >>> query.set_dataset("hsapiens_gene_ensembl")
    >>> query.add_attribute("ensembl_gene_id")
    >>> query.add_attribute("ensembl_peptide_id")
    >>> query.add_filter("chromosome_name", "1")
    >>> count = query.get_count()
"""

import sys, os
import urllib2
import traceback
import posixpath
import xml.dom.minidom as minidom
import shelve

import orngEnviron

class BioMartError(Exception):
    pass

class BioMartServerError(BioMartError):
    pass

class DatabaseNameConflictError(BioMartError):
    pass

class BioMartQueryError(BioMartError):
    pass

def _init__slots_from_line(self, line):
    for attr, val in zip(self.__slots__, line.split("\t")):
        setattr(self, attr, val)
        
class Attribute(object):
    """ A class representing an attribute in a BioMart dataset (returned by BioMartDataset.attributes())
    Attributes::
        - *internalName*    - Internal attribute name
        - *name*    - Human readable name of the attribute
        - *description*    - Human readable description of the attribute
    """
    __slots__ = ["internalName", "name", "description", "_4", "format", "table", "_7"]
    __init__ = _init__slots_from_line
    
    def __repr__(self):
        return "Attribute('%s')" % "\t".join(getattr(self, name, "") for name in self.__slots__)
     
            
class Filter(object):
    """ A class representing a filter for a BioMart dataset (returned by BioMartDataset.filters())
    Attributes::
        - *internalName* - Internal filter name
        - *name*    - Filter name
        - *values*    - Lists possible filter values
        - *description*    - Filter description
    """
    __slots__ =  ["internalName", "name", "values", "description", "type", "qualifiers", "_6", "attribute_table"]
    __init__ = _init__slots_from_line
    
    def __repr__(self):
        return "Filter('%s')" % "\t".join(getattr(self, name, "") for name in self.__slots__)

class xml_node(object):
    def __init__(self, dom):
        self.tag = dom.tagName
        self.attributes = dict(dom.attributes.items())
        self.children = [xml_node(child) for child in dom.childNodes if not isinstance(child, minidom.Text)]
        for child in self.children:
            self.__dict__.setdefault(child.tag, []).append(child)
            
        self.data = "".join([node.data for node in dom.childNodes if isinstance(child, minidom.Text) and node.data.strip()])
            
    def elements(self, tag):
        if self.tag == tag:
            yield self
        for child in self.children:
            for el in child.elements(tag):
                yield el
                
    def __repr__(self):
        return "<%s %s>%s" % (self.tag, self.attributes, self.data) + "\n".join(str(c) for c in self.children)
        
def de_dom(xml, element):
    import xml.dom.minidom as dom
    dom = dom.parse(xml)
    dom.normalize()
    return xml_node(dom.getElementsByTagName(element)[0])

def de_tab(text, sep="\t"):
    return [line.split(sep) for line in text.split("\n") if line.strip()]

def cached(func):
    from functools import wraps
    @wraps(func)
    def f(self):
        if getattr(self, "_cache_" + func.__name__, None) is None:
            setattr(self, "_cache_" + func.__name__, func(self))
        return getattr(self, "_cache_" + func.__name__)
    return f
    
def safe_iter(generator):
    def _iter(generator):
        while True:
            try:
                yield generator.next()
            except StopIteration:
                raise StopIteration
            except Exception, ex:
                print >> sys.stderr, "An error occured during iteration:\n"
                traceback.print_exc(file=sys.stderr)
    return list(_iter(generator))
        
class BioMartConnection(object):
    """ A connection to a BioMart martservice server
    
    Example::
    >>> connection = BioMartConnection("http://www.biomart.org/biomart/martservice")
    >>> registry = connection.registry()
    """
    
    FOLLOW_REDIRECTS = False
    def __init__(self, url=None, cache=None):
        self.url = url if url is not None else "http://www.biomart.org/biomart/martservice"
        self._cache = cache if cache is not None else shelve.open(os.path.join(orngEnviron.bufferDir, "BioMartCache.pck"))
        
    def request_url(self, **kwargs):
        url = self.url + "?" + "&".join("%s=%s" % item for item in kwargs.items() if item[0] != "POST")
#        print >> sys.stderr, url
        return url.replace(" ", "%20")
    
    def request(self, **kwargs):
        """ Request using RESTful access. Use keywords to pass arguments
        Example::
        >>> connection.request(type="registry").read()  
        """
        url = self.request_url(**kwargs)
        if str(url) not in self._cache:
            response = urllib2.urlopen(url).read()
            self._cache[str(url)] = response
            if hasattr(self._cache, "sync"):
                self._cache.sync()
        import StringIO
        return StringIO.StringIO(self._cache[str(url)]) #urllib2.urlopen(self.request_url(**kwargs))
    
    def registry(self, **kwargs):
        """ Return a BioMartRegistry instance
        """
        return BioMartRegistry(self)
    
    def datasets(self, **kwargs):
        text = self.request(type="datasets", **kwargs).read()
        if text.strip().startswith("Mart name conflict"):
            raise DatabaseNameConflictError(text)
        elif text.strip().startswith("Problem retrieving datasets"):
            raise BioMartError(text)
        return de_tab(text)
    
    def attributes(self, **kwargs):
        text = self.request(type="attributes", **kwargs).read()
        if text.strip().startswith("Dataset name conflict"):
            raise DatabaseNameConflictError(text)
        return de_tab(text)
    
    def filters(self, **kwargs):
        test = self.request(type="filters", **kwargs)
        if text.strip().startswith("Dataset name conflict"):
            raise DatabaseNameConflictError(text)
        return de_tab(text)
    
    def configuration(self, **kwargs):
        try:
            return BioMartDatasetConfig(de_dom(self.request(type="configuration", **kwargs), "DatasetConfig"))
        except Exception, ex:
            print >> sys.stderr, ex
            if "virtualSchema" not in kwargs:
                return BioMartDatasetConfig(de_dom(self.request(type="configuration", virtualSchema="default", **kwargs), "DatasetConfig"))
            raise
    
class BioMartRegistry(object):
    """ A class representing a BioMart registry 
    
    Example::
    >>> for schema in registry.virtual_schemas():
    ...    print schema.name
    ...
    """
    def __init__(self, file=None):
        self.connection = file if isinstance(file, BioMartConnection) else None
        if self.connection:
            file = self.connection.request(type="registry")
        else:
            file = open(file, "rb") if isinstance(file, basestring) else file
        self.registry = de_dom(file, "MartRegistry")
    
    @cached
    def virtual_schemas(self):
        """ Return a list of BioMartVirtualSchema instances representing each schema
        """
        schemas = [schema.attributes.get(name, "default") for name in self.registry.elements("virtualSchema")]
        if not schemas:
            schemas = [BioMartVirtualSchema(self.registry, name="default", connection=self.connection)]
        else:
            schemas = [BiomartVirtualSchema(schema, name=schema.attributes.get("name", "default"),
                            connection=self.connection) for schema in self.registry.elements("virtualSchema")]
        return schemas
        
    @cached
    def marts(self):
        """ Return a list off all 'mart' instances (BioMartDatabase instances) regardless of their virtual schemas  
        """
        return reduce(list.__add__, safe_iter(schema.marts() for schema in self.virtual_schemas()), [])
    
    def databases(self):
        """ Same as marts()
        """
        return self.marts()
    
    def datasets(self):
        """ Return a list of all datasets (BioMartDataset instances) from all marts regardless of their virtual schemas
        """
        return reduce(list.__add__, safe_iter(mart.datasets() for mart in self.marts()), [])
    
    def query(self, **kwargs):
        """ Return an initialized BioMartQuery object with registry set to self.
        Pass additional arguments to BioMartQuery.__init__ with keyword arguments
        """
        return BioMartQuery(self, **kwargs)
    
    def links_between(self, exporting, importing, virtualSchema="default"):
        """ Return all links between exporting and importing datasets in the virtualSchema
        """
        schema = [schema for schema in self.virtual_schemas() if schema.name == virtualSchema][0]
        return schema.links_between(exporting, importing)
    
    def get_path(self, exporting, importing):
        pass
    
    def __iter__(self):
        return iter(self.marts())
    
    def _find_mart(self, name):
        try:
            return [mart for mart in self.marts() if name in [getattr(mart, "name", None),
                                    getattr(mart, "internalName", None)]][0]
        except IndexError, ex:
            raise ValueError(name)
        
    def __getitem__(self, name):
        return self._find_mart(name)
    
    def search(self, string, relevance=False):
        results = []
        for mart in self.marts():
            results.extend([(rel, mart, dataset) for rel, dataset in  mart.search(string, True)])
        results = sorted(results, reverse=True)
        if not relevance:
            results = [(mart, dataset) for _, mart, dataset in results]
        return results
    
class BioMartVirtualSchema(object):
    """ A class representation of a virtual schema.
    """
    def __init__(self, locations=None, name="default", connection=None):
        self.locations = locations
        self.name = name
        self.connection = connection
        
    @cached
    def marts(self):
        """ Return a list off all 'mart' instances (BioMartDatabase instances) in this schema
        """
        return safe_iter(BioMartDatabase(connection=self.connection, **dict((str(key), val) \
                for key, val in loc.attributes.items())) for loc in self.locations.elements("MartURLLocation"))
    
    def databases(self):
        """ Same as marts()
        """
        return self.marts()
    
    @cached
    def datasets(self):
        """ Return a list of all datasets (BioMartDataset instances) from all marts in this schema
        """
        return reduce(list.__add__, safe_iter(mart.datasets() for mart in self.marts()), [])
    
    def links_between(self, exporting, importing):
        """ Return a list of link names from exporting dataset to importing dataset 
        """
        exporting = self[exporting]
        importing = self[importing]
        exporting = exporting.configuration().exportables()
        importing = importing.configuration().importables()
        exporting = set([(ex.linkName, getattr(ex, "linkVersion", "")) for ex in exporting])
        importing = set([(im.linkName, getattr(ex, "linkVersion", "")) for im in importing])
        links = exporting.intersection(importing)
        return [link[0] for link in links]
    
    def links(self):
        """ Return a list of (linkName, linkVersion) tuples defined by datasets in the schema
        """
        pass
    
    def query(self, **kwargs):
        """ Return an initialized BioMartQuery object with registry and virtualSchema set to self.
        Pass additional arguments to BioMartQuery.__init__ with keyword arguments
        """
        return BioMartQuery(self, virtualSchema=self.name, **kwargs)

    def __iter__(self):
        return iter(self.marts())
    
    def __getitem__(self, key):
        try:
            return self._find_dataset(key)
        except ValueError:
            raise KeyError(key)
        
    def _find_dataset(self, dataset):
        try:
            index = reduce(list.__add__, safe_iter(mart._datasets_index() for mart in self.marts()), [])
            data = [data for data in index if dataset in [data.get("internalName"), data.get("name")]][0]
            return BioMartDataset(connection=self.connection, **data)
        except IndexError:
            raise ValueError(dataset)
    
class BioMartDatabase(object):
    """ An object representing a BioMart 'mart' instance.
    Arguments::
        - *name*   - Name of the mart instance ('ensembl' by default)
        - *virtualSchema*    - Name of the virtualSchema this dataset belongs to ("default" by default)
        - *connection*    - An optional BioMartConnection instance
    """
    host = "www.biomart.org"
    path = "/biomart/martservice"
    def __init__(self, name="ensembl", virtualSchema="default", connection=None, **kwargs):
#        self.con = url if isinstance(url, BioMartConnection) else BioMartConnection(url)
        self.name = name
        self.virtualSchema = virtualSchema
        self.__dict__.update(kwargs.items())
        self.connection = BioMartConnection("http://" + self.host + ":" + getattr(self, "port", "80") + self.path) if connection is None \
                            or (kwargs.get("redirect", None) == "1" and BioMartConnection.FOLLOW_REDIRECTS) else connection

    @cached    
    def _datasets_index(self):
        keys = ["dataset_type", "internalName", "displayName", "visible", "_4", "_5", "_6", "virtualSchema", "date"]
        try:
            datasets = self.connection.datasets(mart=self.name, virtualSchema=self.virtualSchema)
        except BioMartError, ex:
            if self.virtualSchema == "default":
                print >> sys.stderr, "error", ex
                datasets = self.connection.datasets(mart=self.name)
            else:
                raise
        return [dict(zip(keys, line)) for line in datasets]
    
    @cached
    def datasets(self):
        """ Return a list of all datasets (BioMartDataset instances) in this database
        """
        return safe_iter(BioMartDataset(connection=self.connection, **dataset) for dataset in self._datasets_index())
    
    def dataset_attributes(self, dataset, **kwargs):
        """ Return a list of dataset attributes
        """
        dataset = self._find_dataset(dataset)
        return dataset.attributes()
    
    def dataset_filters(self, dataset, **kwargs):
        """ Return a list of dataset filters
        """
        dataset = self._find_dataset(dataset)
        return dataset.filters()
    
    def dataset_query(self, dataset, filters=[], attributes=[]):
        """ Return an dataset query based on dataset, filters and attributes
        """
        dataset = self._find_dataset(dataset)
        return BioMartQuery(self.con, [dataset], filters, attributes).run()
        
    def _find_dataset(self, dataset):
        try:
            data = [data for data in self._datasets_index() if dataset in [data.get("internalName"), data.get("displayName")]][0]
            return BioMartDataset(connection=self.connection, **data)
        except IndexError:
            raise ValueError(dataset)
        
    def __getitem__(self, name):
        return self._find_dataset(name)
    
    def search(self, string, relevance=False):
        strings = string.lower().split()
        results = []
        for dataset, conf in safe_iter((dataset, dataset.configuration()) for dataset in self.datasets()):
            trees = conf.attribute_pages() + conf.attribute_groups() + \
                    conf.attribute_collections() + conf.attributes()
                    
            names = " ".join([getattr(tree, "displayName", "") for tree in trees]).lower()
            count = sum([names.count(s) for s in strings])
            if count:
                results.append((count ,dataset))
        return results
        
class BioMartDataset(object):
    """ An object representing a BioMart dataset (returned by BioMartDatabase)
    """
    def __init__(self, mart="ensembl", name="hsapiens_gene_ensembl", virtualSchema="default", connection=None, **kwargs):
        self.conenction = connection if connection is not None else BioMartConnection()
        self.mart = mart
        self.name = name
        self.virtualSchema = virtualSchema
        self.__dict__.update(kwargs)
        self.internalName = kwargs.get("internalName", name)
        self._attributes = None
        self._filters = None
        
    def attributes(self):
        """ Return a list of available attributes for this dataset (Attribute instances)
        """
        if self._attributes is None:
            attrs = self.conenction.request(type="attributes", dataset=self.internalName,
                                        virtualSchema=self.virtualSchema).read().split("\n")
            self._attributes = [Attribute(line) for line in attrs if line.strip()]
        return self._attributes
    
    def filters(self):
        """ Return a list of available filters for this dataset (Filter instances)
        """
        if self._filters is None:
            filters = self.conenction.request(type="filters", dataset=self.internalName,
                                        virtualSchema=self.virtualSchema).read().split("\n")
            self._filters = [Filter(line) for line in filters if line.strip()]
        return self._filters
    
    def configuration(self):
        """ Return the configuration tree for this dataset (BioMartDatasetConfig instance)
        """
        return self.conenction.configuration(dataset=self.internalName, virtualSchema=self.virtualSchema)
            
    def get_data(self, attributes=[], filters=[], unique=False):
        """ Constructs and runs a BioMartQuery returning its results 
        """ 
        return BioMartQuery(self.conenction, dataset=self, attributes=attributes, filters=filters,
                             uniqueRows=unique, virtualSchema=self.virtualSchema).run()
    
    def count(self, filters=[], unique=False):
        """ Constructs and runs a BioMartQuery to count the number of returned lines
        """
        return BioMartQuery(self.conenction, dataset=self, filters=filters, uniqueRows=unique,
                            virtualSchema=self.virtualSchema).get_count()
                            
    def get_example_table(self, attributes=[], filters=[], unique=False):
        import orange
        data = self.get_data(attributes, filters, unique)
#        print data
        name = lambda attr: str(getattr(attr, "displayName", getattr(attr, "name", attr)))
        domain = orange.Domain([orange.StringVariable(name(attr)) for attr in attributes], False)
        return orange.ExampleTable(domain, [line.strip().split("\t") for line in data.split("\n") if line.strip()])
        

class BioMartQuery(object):
    """ A class for constructing a query to run on a BioMart server
    Example::
    >>> BioMartQuery(connection, dataset="gene", attributes=["gene_name_", "dictybaseid", "gene_chromosome"],
    ... filters=[("chromosome", 1), ("strand", 1)]).run()
    >>> BioMartQuery(connection, dataset="hsapiens_gene_ensembl", attributes=["ensembl_transcript_id",
    ... "chromosome_name"], filters=[("chromosome_name", ["22"])]).get_count()
    >>> #Equivalent to
    ...
    >>> query = BioMartQuery(connection)
    >>> query.set_dataset("hsapiens_gene_ensembl")
    >>> query.add_filter("chromosome_name", "22")
    >>> query.add_attribute("ensembl_transcript_id")
    >>> query.add_attribute("chromosome_name")
    >>> query.get_count()
    """
    class XMLQuery(object):
        XML = """<?xml version="1.0" encoding="UTF-8"?> 
<!DOCTYPE Query> 
<Query virtualSchemaName = "%(virtualSchemaName)s" formatter = "%(formatter)s" header = "%(header)s" uniqueRows = "%(uniqueRows)s" count = "%(count)s" datasetConfigVersion = "%(version)s" >  
%(datasets_xml)s 
</Query>"""
        version = "0.7"
        def __init__(self, query):
            self.query = query
            self._query = query._query
            self.__dict__.update(query.__dict__) 
        
        def get_xml(self, count=None, header=False):
            datasets_xml = "\n".join(self.query_dataset(*query) for query in self._query)
            
            links_xml = self.query_links(self._query)
            datasets_xml += "\n" + links_xml
            
            count = self.count if count is None else count 
            args = dict(datasets_xml=datasets_xml,
                        uniqueRows="1" if self.uniqueRows else "0",
                        count="1" if count else "",
                        header= "1" if header else "0",
                        virtualSchemaName=self.virtualSchema,
                        formatter=self.format,
                        version=self.version)
            xml = self.XML % args
            return xml
        
        def query_attributes(self, attributes):
            def name(attr):
                if isinstance(attr, Attribute):
                    return attr.internalName
                elif isinstance(attr, BioMartAttribute):
                    return attr.internalName
                else:
                    return str(attr)
                
            xml = "\n".join('\t\t<Attribute name = "%s" />' % name(attr) for attr in attributes)
            return xml
        
        def query_filter(self, filters):
            def name(filter):
                if isinstance(filter, Filter):
                    return filter.internalName
                elif isinstance(filter, BioMartFilter):
                    return getattr(filter, "field", getattr(filter, "internalName"))
                else:
                    return str(filter)
                
            value = lambda value: str(value) if not isinstance(value, list) else ",".join([str(v) for v in value])
            xml = "\n".join('\t\t<Filter name = "%s" value="%s" />' % (name(fil), value(val)) for fil, val in filters)
            return xml
        
        def query_dataset(self, dataset, attributes, filters):
            name, interface = (dataset.internalName, "default") if \
                    isinstance(dataset, BioMartDataset) else (dataset, "default")
            xml = '\t<Dataset name = "%s" interface = "%s" >\n' % (name, interface)
            xml += self.query_attributes(attributes) + "\n"
            xml += self.query_filter(filters)
            xml += '\n\t</Dataset>'
            return xml
        
        def query_links(self, query):
            def name(dataset):
                if isinstance(dataset, BioMartDataset):
                    return getattr(dataset, "internalName", getattr(dataset, "name", None))
                else:
                    return str(dataset)
                
            if len(query) == 2:
                return '\t<Links source="%s" target="%s"/>' % (name(query[0][0]), name(query[1][0]))
            else:
                return ""
        
    class XMLQueryOld(XMLQuery):
        version = "0.4"
        def query_filter(self, filters):
            def name(filter):
                if isinstance(filter, Filter):
                    return filter.internalName
                elif isinstance(filter, BioMartFilter):
                    return getattr(filter, "field", getattr(filter, "internalName"))
                else:
                    return str(filter)
                
            value = lambda value: str(value) if not isinstance(value, list) else ",".join([str(v) for v in value])
            xml = "\n".join('\t\t<ValueFilter name = "%s" value="%s" />' % (name(fil), value(val)) for fil, val in filters)
            return xml
        
    def __init__(self, registry, virtualSchema="default", dataset=None, attributes=[], filters=[], count=False, uniqueRows=False,
                 format="TSV"):
        if isinstance(registry, BioMartConnection):
            self.registry = registry.registry()
            self.virtualSchema = virtualSchema
        elif isinstance(registry, BioMartVirtualSchema):
            self.registry = registry.connection.registry()
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
        stream = self.registry.connection.request(query=self.xml_query(count=count, header=header).replace("\n", "").replace("\t", ""))
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
        dataset = dataset.internalName if isinstance(dataset, BioMartDataset) else dataset
        
        conf = self.registry.connection.configuration(dataset=dataset, virtualSchema=self.virtualSchema)
        self.version = version = getattr(conf, "softwareVersion", "0.4")
        
        if version > "0.4":
            xml = self.XMLQuery(self).get_xml(count, header)
        else:
            xml = self.XMLQueryOld(self).get_xml(count, header)
        return xml
    
    def xml_query_attributes(self, attributes):
        def name(attr):
            if isinstance(attr, Attribute):
                return attr.internalName
            elif isinstance(attr, BioMartAttribute):
                return attr.internalName
            else:
                return str(attr)
            
        xml = "\n".join('\t\t<Attribute name = "%s" />' % name(attr) for attr in attributes)
        return xml
    
    def xml_query_filter(self, filters):
        def name(filter):
            if isinstance(filter, Filter):
                return filter.internalName
            elif isinstance(filter, BioMartFilter):
                return getattr(filter, "field", getattr(filter, "internalName"))
            else:
                return str(filter)
            
        value = lambda value: str(value) if not isinstance(value, list) else ",".join([str(v) for v in value])
        xml = "\n".join('\t\t<ValueFilter name = "%s" value="%s" />' % (name(fil), value(val)) for fil, val in filters)
        return xml
    
    def xml_query_dataset(self, dataset, attributes, filters):
        name, interface = (dataset.internalName, "default") if \
                isinstance(dataset, BioMartDataset) else (dataset, "default")
        xml = '\t<Dataset name = "%s" interface = "%s" >\n' % (name, interface)
        xml += self.xml_query_attributes(attributes) + "\n"
        xml += self.xml_query_filter(filters)
        xml += '\n\t</Dataset>'
        return xml
    
    def xml_query_links(self, query):
        def name(dataset):
            if isinstance(dataset, BioMartDataset):
                return getattr(dataset, "internalName", getattr(dataset, "name", None))
            else:
                return str(dataset)
            
        if len(query) == 2:
            return '\t<Links source="%s" target="%s"/>' % (name(query[0][0]), name(query[1][0]))
        else:
            return ""
        
    def get_example_table(self):
        import orange
        data = self.run(count=False, header=True)
        
        if self.format.lower() == "tsv":
            header, data = data.split("\n", 1)
            domain = orange.Domain([orange.StringVariable(name) for name in header.split("\t")], False)
            data = [line.split("\t") for line in data.split("\n") if line.strip()]
            return orange.ExampleTable(domain, data) if data else None
        elif self.format.lower() == "fasta":
            domain = orange.Domain([orange.StringVariable("id"), orange.StringVariable("sequence")], False) #TODO: meaningful id
            examples = []
            from StringIO import StringIO
            from Bio import SeqIO
            for seq in SeqIO.parse(StringIO(data)):
                examples.append((seq.id, str(seq.seq)))
            return orange.ExampleTable(domain, examples)
        else:
            raise BioMartError("Unsupported format: %" % self.format)


class BioMartConfigurationTree(object):
    _configuration = {}
    def __init__(self, tree):
        self._tree = tree
        self.__dict__.update(tree.attributes.items())
        
    def _get_by_tag(self, tag, name=None):
        return [self._factory(tree) for tree in self._tree.elements(tag)]
    
    def _get_by_name(self, tag, name=None):
        trees = [tree for tree in self._get_by_tag(tag) if name in [getattr(tree, "name", None),
                                                                    getattr(tree, "internalName", None)]]
        return trees[0] if trees else None
    
    def _get_children_by_tag(self, tag, name=None):
        trees = getattr(self._tree, tag, [])
#        print tag, trees
        return [self._factory(tree) for tree in trees]
        
    @classmethod
    def _factory(cls, tree):
        return cls._configuration.get(tree.tag)(tree)

_configuration_names = {
"dataset_config": "DatasetConfig",
"filter_page": "FilterPage",
"filter_group": "FilterGroup",
"filter_collection": "FilterCollection",
"filter": "FilterDescription",
"option": "Option",
"push_action": "PushAction",
"attribute_page": "AttributePage",
"attribute_group": "AttributeGroup",
"attribute_collection": "AttributeCollection",
"attribute": "AttributeDescription",
"exportable": "Exportable",
"importable": "Importable",
}

def _configuration(tagName, names = []):
    def configurator(cls):
        cls._configuration[_configuration_names[tagName]] = cls
        for name in names:
            setattr(cls, name + "s", lambda self, name=name: self._get_by_tag(_configuration_names[name]))
            setattr(cls, name, lambda self, search, name=name: self._get_by_name(_configuration_names[name], search))
        return cls
    return configurator

class BioMartDatasetConfig(BioMartConfigurationTree):
    pass
_configuration("dataset_config", [name for name in _configuration_names if name != "dataset_config"])(BioMartDatasetConfig)

class BioMartExportable(BioMartConfigurationTree):
    pass
_configuration("exportable")(BioMartExportable)

class BioMartImportable(BioMartConfigurationTree):
    pass
_configuration("importable")(BioMartImportable)

class BioMartFilterPage(BioMartConfigurationTree):
    pass
_configuration("filter_page", ["filter_group", "filter_collection", "filter"])(BioMartFilterPage)

class BioMartFilterGroup(BioMartConfigurationTree):
    pass
_configuration("filter_group", ["filter_collection", "filter"])(BioMartFilterGroup)

class BioMartFilterCollection(BioMartConfigurationTree):
    pass
_configuration("filter_collection", ["filter"])(BioMartFilterCollection)

class BioMartFilter(BioMartConfigurationTree):
    def is_pointer(self):
        return hasattr(self, "pointerAttribute") or hasattr(self, "pointerFilter")
    def options(self):
        return self._get_children_by_tag("Option")
_configuration("filter", [])(BioMartFilter)

class BioMartOption(BioMartConfigurationTree):
    def push_actions(self):
        return self._get_children_by_tag("PushAction")
    def options(self):
        return self._get_children_by_tag("Option")
_configuration("option", [])(BioMartOption)

class BioMartPushAction(BioMartConfigurationTree):
    def options(self):
        return self._get_children_by_tag("Option")
_configuration("push_action", [])(BioMartOption)

class BioMartAttributePage(BioMartConfigurationTree):
    pass
_configuration("attribute_page", ["attribute_group", "attribute_collection", "attribute"])(BioMartAttributePage)

class BioMartAttributeGroup(BioMartConfigurationTree):
    pass
_configuration("attribute_group", ["attribute_collection", "attribute"])(BioMartAttributeGroup)

class BioMartAttributeCollection(BioMartConfigurationTree):
    pass
_configuration("attribute_collection", ["attribute"])(BioMartAttributeCollection)

class BioMartAttribute(BioMartConfigurationTree):
    def is_pointer(self):
        return hasattr(self, "pointerAttribute") or hasattr(self, "pointerFilter")
    def target(self, connection=None):
        connection = connection if connection else BioMartConnection()
        dataset = self.pointerDataset
        if hasattr(self, "pointerAttribute"):
            name = self.pointerAttribute
            attr = connection.configuration(dataset=dataset)
         
    
_configuration("attribute")(BioMartAttribute)
    

if __name__ == "__main__":
    con = BioMartConnection("http://www.biomart.org/biomart/martservice")
#    for schema in  con.registry().virtual_schemas():
#        for mart in schema.databases()[:1]:
#            for dataset in mart.datasets()[:1]:
#                print dataset.attributes()
 
    database = BioMartDatabase(name="dicty", connection=con)
    datasets = database.datasets()
    print datasets
    dataset = datasets[2]
    configuration =  dataset.configuration()
    attr = dataset.attributes()
    filters = dataset.filters()
    reg = con.registry()
#    print reg.links_between("dictygo", "gene")
    print BioMartQuery(con, dataset=dataset, attributes=attr[:2], format="FASTA").get_count()
    query = BioMartQuery(con, dataset="gene", attributes=["seq_dictybaseid"], filters=[("chromosome", 1), ("strand", 1)], format="FASTA")
    query.set_dataset("dna")
    query.add_attribute("upstream_intergenic_raw")
    query.add_filter("upstream_flank", 100)
    query.add_filter("downstream_flank", 100)
    print query.run()
    
