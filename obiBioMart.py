
import sys
import urllib2
import traceback
import posixpath
import xml.dom.minidom as minidom

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
    __slots__ = ["internal_name", "name", "description", "page_internal_name", "_4", "_5", "interface"]
    __init__ = _init__slots_from_line
    
    def __repr__(self):
        return "\t".join(getattr(self, name, "") for name in self.__slots__)
            
class Filter(object):
    __slots__ =  ["internal_name", "name", "values", "description", "filter_type", "filter_method", "_6", "_7"]
    __init__ = _init__slots_from_line
    
    def __repr__(self):
        return "\t".join(getattr(self, name, "") for name in self.__slots__)

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
#    print xml
    dom = dom.parse(xml)
    dom.normalize()
    return xml_node(dom.getElementsByTagName(element)[0])

def de_tab(text, sep="\t"):
    return [line.split(sep) for line in text.split("\n") if line.strip()]

def cached(func):
    from functools import wraps
    @wraps
    def f(self):
        if getattr(self, "_" + func.__name__, None) is None:
            setattr(self, "_" + func.__name__, func(self))
        return getattr(self, "_" + func.__name__)
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
                traceback.print_last(file=sys.stderr)
    return list(_iter(generator))
    
        
class BioMartConnection(object):
    """ A connection to a martservice server
    """
    def __init__(self, url=None):
        self.url = url if url is not None else "http://www.biomart.org/biomart/martservice"
        
    def request(self, **kwargs):
        url = self.url + "?" + "&".join("%s=%s" % item for item in kwargs.items() if item[0] != "POST")
        print url
        url = url.replace(" ", "%20")
        return urllib2.urlopen(url) #kwargs.get("POST", None))
    
    def registry(self, **kwargs):
        return BioMartRegistry(self)
    
    def datasets(self, **kwargs):
        text = self.request(type="datasets", **kwargs).read()
        if text.strip().startswith("Mart name conflict"):
            raise DatabaseNameConflictError(text)
        elif text.strip().startswith("Problem retrieving datasets"):
            raise BioMartError(text)
        return de_tab(self.request(type="datasets", **kwargs).read())
    
    def attributes(self, **kwargs):
        text = self.request(type="attributes", **kwargs).read()
        if text.strip().startswith("Dataset name conflict"):
            raise DatabaseNameConflictError(text)
        return de_tab(text)
    
    def filters(self, **kwargs):
        return self.request(type="filters", **kwargs)
    
    def configuration(self, **kwargs):
        return BioMartDatasetConfig(de_dom(self.request(type="configuration", **kwargs), "DatasetConfig"))
    
class BioMartRegistry(object):
    def __init__(self, file=None):
        self.connection = file if isinstance(file, BioMartConnection) else None
        if self.connection:
            file = self.connection.request(type="registry")
        else:
            file = open(file, "rb") if isinstance(file, basestring) else file
        self.registry = de_dom(file, "MartRegistry")
    
#    @cached
    def virtual_schemas(self):
        schemas = [schema.attributes.get(name, "default") for name in self.registry.elements("virtualSchema")]
        if not schemas:
            schemas = [BioMartVirtualSchema(self.registry, name="default", connection=self.connection)]
        else:
            schemas = [BiomartVirtualSchema(schema, name=schema.attributes.get("name", "default"), connection=self.connection) for schema in self.registry.elements("virtualSchema")]
        return schemas
        
    def databases(self):
        """ Same as marts()
        """
        return self.marts()
    
#    @cached
    def marts(self):
        """ Return a list off all 'mart' instances
        """
        return reduce(list.__add__, safe_iter(schema.marts() for schema in self.virtual_schemas()), [])
    
#    @cached
    def datasets(self): #, virtualSchema="default"):
        """ Return all datasets.
        """
        return reduce(list.__add__, safe_iter(mart.datasets() for mart in self.marts()), [])
    
    def links_between(self, exporting, importing, virtualSchema="default"):
        """ Return all links between exporting and importing datasets in the virtualSchema
        """
        schema = [schema for schema in self.virtual_schemas() if schema.name == virtualSchema][0]
        return schema.links_between(exporting, importing)
    
    def __iter__(self):
        return iter(self.marts())
    
class BioMartVirtualSchema(object):
    def __init__(self, locations=None, name="default", connection=None):
        self.locations = locations
        self.name = name
        self.connection = connection
        
    def marts(self):
        """ Return all 'mart' instances belonging to this shema
        """
        return safe_iter(BioMartDatabase(connection=self.connection, **dict((str(key), val) for key, val in loc.attributes.items())) for loc in self.locations.elements("MartURLLocation"))
    
    def databases(self):
        """ Same as marts()
        """
        return self.marts()
    
    def datasets(self):
        """ Return all datasets in this virtual schema
        """
        return reduce(list.__add__, safe_iter(mart.datasets() for mart in self.marts()), [])
    
    def links_between(self, exporting, importing):
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

    def __iter__(self):
        return iter(self.marts())
    
    def __getitem__(self, key):
        try:
            return self._find_dataset(key)
        except ValueError:
            raise KeyError(key)
        
    def _find_dataset(self, dataset):
        try:
            return [data for data in self.datasets() if dataset in [data, data.internal_name, data.name]][0]
        except IndexError:
            raise ValueError(dataset)
    
class BioMartDatabase(object):
    """ An object representing a BioMart 'mart' instance.
    Arguments::
        - *name*   name of the mart instance ('biomart' by default)
        - *virtualSchema*    name of the virtualSchema this dataset belongs to ("default" by default) 
    """
    def __init__(self, name="biomart", virtualSchema="default", connection=None, **kwargs):
#        self.con = url if isinstance(url, BioMartConnection) else BioMartConnection(url)
        self.name = name
        self.virtualSchema = virtualSchema
        self.__dict__.update(kwargs.items())
        self._datasets = None
        self.connection = BioMartConnection("http://" + self.host + ":" + getattr(self, "port", "80") + self.path) if connection is None \
                            or kwargs.get("redirect", None) == "1" else connection
    
    def datasets(self):
        """ Return a list of all datasets in a BioMart database
        """
        if self._datasets is None:
            keys = ["dataset_type", "internal_name", "name", "dataset_visible", "_4", "_5", "_6", "virtual_schema_name", "date"]
#            datasets = safe_iter(self.connection.datasets(mart=self.name)
            datasets = [dict(zip(keys, line)) for line in self.connection.datasets(mart=self.name, virtualSchema=self.virtualSchema)]
            self._datasets = safe_iter(BioMartDataset(dataset, connection=self.connection) for dataset in datasets)
        return self._datasets
    
    
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
            return [data for data in self.datasets() if dataset in [data, data.internal_name, data.name]][0]
        except IndexError:
            raise ValueError(dataset)
        
    def __getitem__(self, name):
        return self._find_dataset(dataset)
        
class BioMartDataset(object):
    """ An object representing a BioMart dataset (returned by BioMartDatabase)
    """
    def __init__(self, dataset_info, connection=None):
        self.conenction = connection
        self.__dict__.update(dataset_info.items())
        self._attributes = None
        self._filters = None
        
    def attributes(self):
        """ Return a list of available attributes for this dataset
        """
        if self._attributes is None:
            attrs = self.conenction.request(type="attributes", dataset=self.internal_name,
                                        virtualSchema=self.virtual_schema_name).read().split("\n")
            self._attributes = [Attribute(line) for line in attrs if line.strip()]
        return self._attributes
    
    def filters(self):
        """ Return a list of available filters for this dataset
        """
        if self._filters is None:
            filters = self.conenction.request(type="filters", dataset=self.internal_name,
                                        virtualSchema=self.virtual_schema_name).read().split("\n")
            self._filters =  [Filter(line) for line in filters if line.strip()]
        return self._filters
    
    def configuration(self):
        """ Return the configuration tree for this dataset 
        """
        return self.conenction.configuration(dataset=self.internal_name)
            
    def get_data(self, attributes=[], filters=[], unique=False):
        """ Constructs and runs a query returning its results 
        """ 
        return BioMartQuery(self.conenction, [self], attributes, filters, unique=unique).run()
    
    def count(self, filters=[], unique=False):
        """ Constructs and runs a query to count the number of returned lines
        """
        return BioMartQuery(self.conenction, [self], filters, unique=unique).get_count()

class BioMartQuery(object):
    """ A class for constructing a query to run on a BioMart server
    Example::
    >>> BioMartQuery(con, ["gene"], attributes=["gene_name_", "dictybaseid", "gene_chromosome"], filters=[("chromosome", 1), ("strand", 1)]).run()
    >>> BioMartQuery(con, ["hsapiens_gene_ensembl"], attributes=["ensembl_transcript_id", "chromosome_name"], filters=[("chromosome_name", ["22"])]).get_count()
    """
    XML = """<?xml version="1.0" encoding="UTF-8"?> 
<!DOCTYPE Query> 
<Query virtualSchemaName = "%(virtual_schema_name)s" formatter = "%(formatter)s" header = "0" uniqueRows = "%(unique_rows)s" count = "%(count)s" datasetConfigVersion = "0.4" > 
%(datasets_xml)s 
</Query>""" 
    def __init__(self, connection, datasets=[], attributes=[], filters=[], count=False, unique=False, format="TSV", virtual_schema_name="default"):
        self.connection = connection
        self.datasets = datasets
        self.attributes = attributes
        self.filters = filters
        self.count = count
        self.unique = unique
        self.format = format
        self.virtual_schema_name = virtual_schema_name
        
    def set_dataset(self, dataset):
        self.datasets.append(dataset)
    
    def add_filter(self, filter, value):
        self.filters.append((filter, value))
    
    def add_attribute(self, attribute):
        self.attributes.append(attribute)
    
    def get_count(self):
        tmp, self.count = self.count, True
        count = self.run()
        count = int(count.strip())
        self.count = tmp
        return count
    
    def run(self):
        stream = self.connection.request(query=self.xml_query().replace("\n", "").replace("\t", ""))
        stream = stream.read()
        if stream.startswith("Query ERROR:"):
            raise BioMartQueryError(stream)
        return stream
    
    def set_unique_rows(self, unique=False):
        self.unique = unique
    
    def xml_query(self):
        datasets_xml = self.xml_query_dataset(self.datasets[0], self.attributes, self.filters)
        args = dict(datasets_xml=datasets_xml,
                    unique_rows="1" if self.unique else "0",
                    count="1" if self.count else "",
                    virtual_schema_name=self.virtual_schema_name,
                    formatter=self.format)
        xml = self.XML % args
        return xml
    
    def xml_query_attributes(self, attributes):
        name = lambda attr: attr.internal_name if isinstance(attr, Attribute) else str(attr)
        xml = "\n".join('\t\t<Attribute name = "%s" />' % name(attr) for attr in attributes)
        return xml
    
    def xml_query_filter(self, filters):
        name = lambda filter: filter.internal_name if isinstance(filter, Filter) else str(filter)
        value = lambda value: str(value) if not isinstance(value, list) else ",".join([str(v) for v in value])
        xml = "\n".join('\t\t<ValueFilter name = "%s" value="%s" />' % (name(fil), value(val)) for fil, val in filters)
        return xml
    
    def xml_query_dataset(self, dataset, attributes, filters):
        name, interface = (dataset.internal_name, "default") if \
                isinstance(dataset, BioMartDataset) else (dataset, "default")
        xml = '\t<Dataset name = "%s" interface = "%s" >\n' % (name, interface)
        xml += self.xml_query_attributes(attributes) + "\n"
        xml += self.xml_query_filter(filters)
        xml += '\n\t</Dataset>'
        return xml
 
 
class BioMartConfigurationTree(object):
    _configuration = {}
    def __init__(self, tree):
        self._tree = tree
        self.__dict__.update(tree.attributes.items())
        
    def _get_by_tag(self, tag, name=None):
        return [self._factory(tree) for tree in self._tree.elements(tag)]
    
    def _get_by_name(self, tag, name=None):
        return [tree for tree in self._get_by_tag(tag) if name in [getattr(tree, "name", None),
                                                                   getattr(tree, "internalName", None)]]
        
    @classmethod
    def _factory(cls, tree):
        print tree.tag
        return cls._configuration.get(tree.tag)(tree)

_configuration_names = {
"dataset_config": "DatasetConfig",
"filter_page": "FilterPage",
"filter_group": "FilterGroup",
"filter_collection": "FilterCollection",
"filter": "FilterDescription",
"option": "Option",
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
        print cls, names
        for name in names:
            setattr(cls, name + "s", lambda self, name=name: self._get_by_tag(_configuration_names[name]))
            setattr(cls, name, lambda self, search, name=name: self._get_by_name(_configuration_names[name], search))
#            print getattr(cls, name)
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
    pass
_configuration("filter", ["option"])(BioMartFilter)

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
    pass
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
    print BioMartQuery(con, [dataset], attr[:2], count=True).get_count().read()
    print BioMartQuery(con, [dataset], ["ensembl_transcript_id"], filters=[("with_go_molecular_function", "")]).get_count().read()
    
