""" 
Python module for manipulating OBO Ontology files (http://www.obofoundry.org/)

TODO:
    - handle escape characters !!!

Example::
    >>> term = OBOObject("Term", id="foo:bar", name="Foo bar")
    >>> print term
    [Term]
    id: foo:bar
    name: Foo bar
    
    >>> ontology = OBOOntology()
    >>> ontology.addObject(term)
    
"""

from itertools import chain
from collections import defaultdict
import itertools

builtinOBOObjects = [
"""[Typedef]
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
definition: Indicates that a term is the intersection of several others [OBO:defs]"""
]

def _splitAndStrip(string, sep):
    head, tail = string.split(sep, 1)
    return head.rstrip(" "), tail.lstrip(" ")


class OBOObject(object):
    """ Represents a generic OBO object (e.g. Term, Typedef, Instance, ...)
    Example::
        >>> term = OBOObject(stanzaType="Term", id="FOO:001", name="bar")
    """
    def __init__(self, stanzaType="Term", **kwargs):
        """ Init from keyword arguments.
        Example::
            >>> term = OBOObject(stanzaType="Term", id="FOO:001", name="bar", def_="Example definition {modifier=frob} ! Comment")
            >>> term = OBOObject(stanzaType="Term", id="FOO:001", name="bar", def_=("Example definition", [("modifier", "frob")], "Comment"))
            >>> term = OBOObject(stanzaType="Term", id="FOO:001", name="bar", def_=("Example definition", [("modifier", "frob")])) # without the comment
            >>> term = OBOObject(stanzaType="Term", id="FOO:001", name="bar", def_=("Example definition",)) # without the modifiers and comment
        """
        self.stanzaType = stanzaType
        
        self.modifiers = []
        self.comments = []
        self.tagValues = []
        self.values = {}
        
        sortedTags = sorted(kwargs.iteritems(), key=lambda key_val: chr(1) if key_val[0] == "id" else key_val[0])
        for tag, value in sortedTags:
            if isinstance(value, basestring):
                tag, value, modifiers, comment = self.parseTagValue(self.name_demangle(tag), value)
            elif isinstance(value, tuple):
                tag, value, modifiers, comment = ((self.name_demangle(tag),) + value + (None, None))[:4]
            self.addTag(tag, value, modifiers, comment)
        
        self.related = set()
        self.relatedTo = set()
            
    @property
    def is_annonymous(self):
        value = self.getValue("is_annonymous")
        return bool(value)
    
    def name_mangle(self, tag):
        """ Mangle tag name if it conflicts with python keyword
        Example::
            >>> term.name_mangle("def"), term.name_mangle("class")
            ('def_', 'class_')
        """
        if tag in ["def", "class", "in", "not"]:
            return tag + "_"
        else:
            return tag
        
    def name_demangle(self, tag):
        """ Reverse of name_mangle
        """
        if tag in ["def_", "class_", "in_", "not_"]:
            return tag[:-1]
        else:
            return tag
        
    def addTag(self, tag, value, modifiers=None, comment=None):
        """ Add `tag`, `value` pair to the object with optional modifiers and
        comment.
        Example::
            >>> term = OBOObject("Term")
            >>> term.addTag("id", "FOO:002", comment="This is an id")
            >>> print term
            [Term]
            id: FOO:002 ! This is an id
             
        """
        self.tagValues.append((tag, value))
        self.modifiers.append(modifiers)
        self.comments.append(comment)
        self.values.setdefault(tag, []).append(value)
        
        #  TODO: fix multiple tags grouping
        if hasattr(self, tag):
            if isinstance(getattr(self, tag), list):
                getattr(self, tag).append(value)
            else:
                setattr(self, tag, [getattr(self, tag)] + [value])
        else:
            setattr(self, self.name_mangle(tag), value)
            
    def update(self, other):
        """ Update the term with tag value pairs from `other` 
        (a OBOObject instance). The tag value pairs are appended
        to the end except for the `id` tag.
        """ 
        for (tag, value), modifiers, comment in zip(other.tagValues, other.modifiers, other.comments):
            if tag != "id":
                self.addTag(tag, value, modifiers, comment)
        
    def getValue(self, tag, group=True):
        if group:
            pairs = [pair for pair in self.tagValues if pair[0] == tag]
            return pairs
        else:
            tag = self.name_mangle(tag)
            if tag in self.__dict__:
                return self.__dict__[tag]
            else:
                raise ValueError("No value for tag: %s" % tag)
        
    def tagCount(self):
        """ Retrun the number of tags in this object
        """
        return len(self.tagValues)
    
    def tags(self):
        """ Retrun an iterator over the (tag, value) pairs.
        """
        for i in range(self.tagCount()):
            yield self.tagValues[i] + (self.modifiers[i], self.comments[i])
        
    def formatSingleTag(self, index):
        """Return a formated string representing index-th tag pair value
        Example::
            >>> term = OBOObject("Term", id="FOO:001", name="bar", def_="Example definition {modifier=frob} ! Comment")
            >>> term.formatSingleTag(0)
            'id: FOO:001'
            >>> term.formatSingleTag(1)
            'def: Example definition { modifier=frob } ! Comment'
        """
        tag, value = self.tagValues[index]
        modifiers = self.modifiers[index]
        comment = self.comments[index]
        res = ["%s: %s" % (tag, value)]
        if modifiers:
            res.append("{ %s }" % modifiers)
        if comment:
            res.append("! " + comment)
        return " ".join(res)
    
    def formatStanza(self):
        """ Return a string stanza representation of this object 
        """
        stanza = ["[%s]" % self.stanzaType]
        for i in range(self.tagCount()):
            stanza.append(self.formatSingleTag(i))
        return "\n".join(stanza)
            
    @classmethod     
    def parseStanza(cls, stanza):
        r""" Parse and return an OBOObject instance from a single stanza.
        Example::
            >>> term = OBOObject.parseStanza("[Term]\nid: FOO:001\nname:bar")
            >>> print term.id, term.name
            FOO:001 bar
            
        """
        lines = stanza.splitlines()
        stanzaType = lines[0].strip("[]")
        tag_values = []
        for line in lines[1:]:
            if ":" in line:
                tag_values.append(cls.parseTagValue(line))
        
        obo = OBOObject(stanzaType)
        for i, (tag, value, modifiers, comment) in enumerate(tag_values):
#            print tag, value, modifiers, comment
            obo.addTag(tag, value, modifiers, comment)
        return obo
    
        
    @classmethod
    def parseTagValue(cls, tagValuePair, *args):
        """ Parse and return a four-tuple containing a tag, value, a list of modifier pairs, comment.
        If no modifiers or comments are present the corresponding entries will be None.
        
        Example::
            >>> OBOObject.parseTagValue("foo: bar {modifier=frob} ! Comment")
            ('foo', 'bar', 'modifier=frob', 'Comment')
            >>> OBOObject.parseTagValue("foo: bar")
            ('foo', 'bar', None, None)
            >>> #  Can also pass tag, value pair already split   
            >>> OBOObject.parseTagValue("foo", "bar {modifier=frob} ! Comment")
            ('foo', 'bar', 'modifier=frob', 'Comment')
        """
        if args and ":" not in tagValuePair:
            tag, rest = tagValuePair, args[0]
        else:
            tag, rest = _splitAndStrip(tagValuePair, ":")
        value, modifiers, comment = None, None, None
        
        if "{" in rest:
            value, rest = _splitAndStrip(rest, "{",)
            modifiers, rest = _splitAndStrip(rest, "}")
        if "!" in rest:
            if value is None:
                value, comment = _splitAndStrip(rest, "!")
            else:
                _, comment = _splitAndStrip(rest, "!")
        if value is None:
            value = rest
            
        if modifiers is not None:
            modifiers = modifiers #TODO: split modifiers in a list
            
        return tag, value, modifiers, comment
         
    def GetRelatedObjects(self):
        """ Obsolete. Use `relatedObjects()` instead
        """
        return self.relatedObjects()
    
    def relatedObjects(self):
        """ Return a list of tuple pairs where the first element is relationship (typedef id)
        is and the second object id whom the relationship applies to.
        """
        result = [(typeId, id) for typeId in ["is_a"] for id in self.values.get(typeId, [])] ##TODO add other defined Typedef ids
        result = result + [tuple(r.split(None, 1)) for r in self.values.get("relationship", [])]
        return result

    def __repr__(self):
        """ Return a string representation of the object in OBO format
        """
        return self.formatStanza()

    def __iter__(self):
        """ Iterates over sub terms
        """
        for typeId, id in self.relatedObjects():
            yield (typeId, id)
        
class Term(OBOObject):
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Term", *args, **kwargs)

class Typedef(OBOObject):
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Typedef", *args, **kwargs)

class Instance(OBOObject):
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Instance", *args, **kwargs)

import re

class OBOOntology(object):
    _RE_TERM = re.compile(r"\[.+?\].*?\n\n", re.DOTALL)
    _RE_HEADER = re.compile(r"^[^[].*?\n\[", re.DOTALL)
    BUILTINS = builtinOBOObjects
    
    def __init__(self):
        """ Init an empty Ontology.
        
        .. note:: Use parseOBOFile to load from a file
        
        """
        self.objects = []
        self.headerTags = []
        self.id2Term = {}
        
    def addObject(self, object):
        """ Add OBOObject instance to  this ontology.
        """
        if object.id in self.id2Term:
            raise ValueError("OBOObject with id: %s already in the ontology" % object.id)
        self.objects.append(object)
        self.id2Term[object.id] = object
        
    def addHeaderTag(self, tag, value):
        """ Add header tag, value pair to this ontology
        """
        self.headerTags.append((tag, value))
        
#    @classmethod
#    def parseOBOFile(cls, file):
#        """ Parse the .obo file and return an OBOOntology instance
#        Example::
#            >>> OBOOntology.parseOBOFile(open("dictyostelium_anatomy.obo", "rb"))
#            <obiOntology.OBOOntology object at ...>
#        """ 
#        ontology = OBOOntology()
#        data = file.read()
#        
#        header = data[:data.index("\n[")]
#        for line in header.splitlines():
#            if line.strip():
#                ontology.addHeaderTag(*line.split(":", 1))
#        
#        imports = [value for  tag, value in ontology.headerTags if tag == "import"]
#        
#        terms = cls.BUILTINS + cls._RE_TERM.findall(data)
#        for term in terms:
#            term = OBOObject.parseStanza(term)
#            ontology.addObject(term)
#            
#        while imports:
#            url = imports.pop(0)
#            imported = self.parseOBOFile(open(url, "rb"))
#            ontology.update(imported)
#        return ontology
    
    @classmethod
    def parseOBOFile(cls, file):
        """ Parse the .obo file and return an OBOOntology instance
        Example::
            >>> OBOOntology.parseOBOFile(open("dictyostelium_anatomy.obo", "rb"))
            <obiOntology.OBOOntology object at ...>
        """ 
        ontology = OBOOntology()
        data = file.read()
        header = data[: data.index("\n[")]
        body = data[data.index("\n[") + 1:]
        for line in header.splitlines():
            if line.strip():
                ontology.addHeaderTag(*line.split(":", 1))
                
        current = None
        #  For speed make these functions local
        startswith = str.startswith
        endswith = str.endswith
        parseTagValue = OBOObject.parseTagValue
        
        builtins = "\n\n".join(cls.BUILTINS)
        for line in itertools.chain(builtins.splitlines(), body.splitlines()):
#            line = line.strip()
            if startswith(line, "[") and endswith(line, "]"):
                current = OBOObject(line.strip("[]"))
            elif startswith(line, "!"):
                pass #  comment
            elif line:
                current.addTag(*parseTagValue(line))
            else: #  empty line is the end of a term
                ontology.addObject(current)
        if current.id not in ontology:
            ontology.addObject(current)
        imports = [value for  tag, value in ontology.headerTags if tag == "import"]
        
        while imports:
            url = imports.pop(0)
            imported = self.parseOBOFile(open(url, "rb"))
            ontology.update(imported)
        return ontology
        
    
    def update(self, other):
        """ Update this ontology with the terms from another. 
        """
        for term in other:
            if term.id in self:
                if not term.is_annonymous:
                    self.term(term.id).update(term)
                else: #  Do nothing
                    pass 
            else:
                self.addObject(term)
                
    def _postLoadProcess(self):
        for obj in self.objects:
            pass
    
    def term(self, id):
        """ Return the OBOObject associated with this id.
        """
        if id in self.id2Term:
            return self.id2Term[id]
        else:
            raise ValueError("Unknown term id: %s" % id)
        
    def terms(self):
        """ Return all `Term` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanzaType == "Term"]
    
    def typedefs(self):
        """ Return all `Typedef` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanzaType == "Typedef"]
    
    def instances(self):
        """ Return all `Instance` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanzaType == "Instance"]
        
    def relatedTerms(self, term):
        """ Return a list of (rel_type, term_id) tuples where rel_type is
        relationship type (e.g. 'is_a', 'has_part', ...) and term_id is the
        id of the term in the relationship.
        """
        term = self.term(term) if not isinstance(term, OBOObject) else term
        related = [(tag, value) for tag in ["is_a"] for value in term.values.get(tag, [])] #TODO: add other typedef ids
        relationships = term.values.get("relationship", [])
        for rel in relationships:
            related.append(tuple(rel.split(None, 1)))
        return related
            
    def toNetwork(self):
        """ Return a orngNetwork instance constructed from this ontology
        """
        edgeTypes = self.edgeTypes()
        terms = self.terms()
        import orngNetwork, orange
        
        network = orngNetwork.Network(len(terms), True, len(edgeTypes))
        network.objects = dict([(term.id, i) for i, term in enumerate(terms)])
        
        edges = defaultdict(set)
        for term in self.terms():
#            related = term.relatedTerms()
            related = self.relatedTerms(term)
            for relType, relTerm in related:
                edges[(term.id, relTerm)].add(relType)
                
        edgeitems = edges.items()
        for (src, dst), eTypes in edgeitems:
            network[src, dst] = [1 if e in eTypes else 0 for e in edgeTypes]
            
        domain = orange.Domain([orange.StringVariable("id"),
                                orange.StringVariable("name"),
                                orange.StringVariable("def"),
                                ], False)
        
        items = orange.ExampleTable(domain)
        for term in terms:
            ex = orange.Example(domain, [term.id, term.name, term.values.get("def", [""])[0]])
            items.append(ex)
        
        relationships = set([", ".join(sorted(eTypes)) for (_, _), eTypes in edgeitems])
        domain = orange.Domain([orange.FloatVariable("u"),
                                orange.FloatVariable("v"),
                                orange.EnumVariable("relationship", values=list(edgeTypes))
                                ], False)
        
        id2index = dict([(term.id, i + 1) for i, term in enumerate(terms)])
        links = orange.ExampleTable(domain)
        for (src, dst), eTypes in edgeitems:
            ex = orange.Example(domain, [id2index[src], id2index[dst], eTypes.pop()])
            links.append(ex)
            
        network.items = items
        network.links = links
        network.optimization = None
        return network
        
    def edgeTypes(self):
        """ Return a list of all edge types in the ontology
        """
        return [obj.id for obj in self.objects if obj.stanzaType == "Typedef"]
        
    def extractSuperGraph(self, terms):
        """ Return all super terms of terms up to the most general one.
        """
        terms = [terms] if type(terms) == str else terms
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(set(id for typeId, id in self[term].related) - visited)
        return visited
    
    def __len__(self):
        return len(self.objects)
    
    def __iter__(self):
        return iter(self.objects)
    
    def __contains__(self, obj):
        return obj in self.id2Term
    
def foundryOntologies():
    """ List ontologies available from the OBOFoundry website
    (`http://www.obofoundry.org/`_) 
    Example::
        >>> foundryOntologies()
        [('Biological process', 'http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/genomic-proteomic/gene_ontology_edit.obo'), ...
    
    """
    import urllib2, re
    stream = urllib2.urlopen("http://www.obofoundry.org/")
    text = stream.read()
    pattern = r'<td class=".+?">\s*<a href=".+?">(.+?)</a>\s*</td>\s*<td class=".+?">.*?</td>\s*<td class=".+?">.*?</td>\s*?<td class=".+?">\s*<a href="(.+?obo)">.+?</a>'
    return re.findall(pattern, text)
    
    
if __name__ == "__main__":
    import doctest
    stanza = '''[Term]
id: FOO:001
name: bar
    '''
    term = OBOObject.parseStanza(stanza)
    doctest.testmod(extraglobs={"stanza": stanza, "term": term}, optionflags=doctest.ELLIPSIS)
        