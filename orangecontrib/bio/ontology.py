"""
==============================
OBO Ontology (:mod:`ontology`)
==============================

This module provides an interface for parsing, creating and manipulating of
`OBO ontologies <http://www.obofoundry.org/>`_.

Construct an ontology from scratch with custom terms ::

    >>> term = OBOObject("Term", id="foo:bar", name="Foo bar")
    >>> print term
    [Term]
    id: foo:bar
    name: Foo bar

    >>> ontology = OBOOntology()
    >>> ontology.add_object(term)
    >>> ontology.add_header_tag("created-by", "ales") # add a header tag
    >>> from StringIO import StringIO
    >>> buffer = StringIO()
    >>> ontology.write(buffer) # Save the ontology to a file like object
    >>> print buffer.getvalue() # Print the contents of the buffer
    created-by: ales
    <BLANKLINE>
    [Term]
    id: foo:bar
    name: Foo bar
    <BLANKLINE>

To load an ontology from a file, pass the file or filename to the
:class:`OBOOntology` constructor or call its load method ::

    >>> buffer.seek(0) # rewind
    >>> ontology = OBOOntology(buffer)
    >>> # Or equivalently
    >>> buffer.seek(0) # rewind
    >>> ontology = OBOOntology()
    >>> ontology.load(buffer)


See the definition of the `.obo file format <http://www.geneontology.org/GO.format.obo-1_2.shtml>`_.

"""

import re
import urllib2
import warnings
import keyword
from collections import defaultdict
from StringIO import StringIO


#: These are builtin OBO objects present in any ontology by default.
BUILTIN_OBO_OBJECTS = [
"""[Typedef]
id: is_a
name: is_a
range: OBO:TERM_OR_TYPE
domain: OBO:TERM_OR_TYPE
definition: The basic subclassing relationship [OBO:defs]""",

"""[Typedef]
id: disjoint_from
name: disjoint_from
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that two classes are disjoint [OBO:defs]""",

"""[Typedef]
id: instance_of
name: instance_of
range: OBO:TERM
domain: OBO:INSTANCE
definition: Indicates the type of an instance [OBO:defs]""",

"""[Typedef]
id: inverse_of
name: inverse_of
range: OBO:TYPE
domain: OBO:TYPE
definition: Indicates that one relationship type is the inverse of another [OBO:defs]""",

"""[Typedef]
id: union_of
name: union_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the union of several others [OBO:defs]""",

"""[Typedef]
id: intersection_of
name: intersection_of
range: OBO:TERM
domain: OBO:TERM
definition: Indicates that a term is the intersection of several others [OBO:defs]"""
]


def _split_and_strip(string, sep):
    """
    Split the `string` by separator `sep` in to two parts and strip
    any whitespace between the inner parts.

    """
    head, tail = _split_esc(string, sep)
    return head.rstrip(" "), tail.lstrip(" ")


def _rsplit_and_strip(string, sep):
    """
    Right split the `string` by separator `sep` in to two parts and
    strip any whitespace between the inner parts.

    """
    head, tail = _rsplit_esc(string, sep)
    return head.rstrip(" "), tail.lstrip(" ")


def _find_esc(string, char):
    i = string.find(char)
    while i != -1:
        if (i > 0 and string[i - 1] != "\\") or string[i - 1] != "\\":
            return i
        else:
            i = string.find(char, i + 1)
    return i


def _rfind_esc(string, char):
    i = string.rfind(char)
    while i != -1:
        if (i > 0 and string[i - 1] != "\\") or string[i - 1] != "\\":
            return i
        else:
            i = string.rfind(char, 0, i - 1)
    return i


def _split_esc(string, sep, _find_esc=_find_esc):
    i = _find_esc(string, sep)
    if i != -1:
        return string[:i], string[i + 1:]
    else:
        return string, ""


def _rsplit_esc(string, sep):
    i = _rfind_esc(string, sep)
    if i != -1:
        return string[:i], string[i + 1:]
    else:
        return string, ""


def parse_tag_value(tag_value_string):
    """
    Parse a tag value string and return a four-tuple containing
    a (tag, value, modifiers, comment). If comment or modifiers are
    not present the corresponding entry will be ``None``.

    >>> parse_tag_value("foo: bar {modifier=frob} ! Comment")
    ('foo', 'bar', 'modifier=frob', 'Comment')
    >>> parse_tag_value("foo: bar")
    ('foo', 'bar', None, None)
    >>> parse_tag_value("foo: bar [baz:0] { fizz=buzz } ! Comment")
    ('foo', 'bar [baz:0]', 'fizz=buzz', 'Comment')

    """
    comment = modifiers = None
    # First get rid of the comment if present
    if _rfind_esc(tag_value_string, "!") != -1:
        tag_value_string, comment = _rsplit_and_strip(tag_value_string, "!")

    # Split on the first unescaped ":"
    tag, value = _split_and_strip(tag_value_string, ":")

    # Split the value on { to get the modifiers if present
    value = value.rstrip()
    if value.endswith("}") and not value.endswith(r"\}") and \
            _rfind_esc(value, "{") != -1:
        value, modifiers = _rsplit_and_strip(value, "{")
        # remove closing } and any whitespace
        modifiers = modifiers[: -1].rstrip()

    return tag, value, modifiers, comment


class OBOObject(object):
    """
    A generic OBO object (e.g. Term, Typedef, Instance, ...).
    Example::

        >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar")

        >>> term = OBOObject(
        ...     stanza_type="Term",
        ...     id="FOO:001",
        ...     name="bar",
        ...     def_="Example definition { modifier=frob } ! Comment"
        ... )
        ...

    An alternative way to specify tags in the constructor::

        >>> term = OBOObject(stanza_type="Term", id="FOO:001", name="bar",
        ...                  def_=("Example definition",
        ...                        [("modifier", "frob")],
        ...                        "Comment"))
        ...

    .. note::
        Note the use of ``def_`` to define the 'def' tag. This is to
        avoid the name clash with the python's ``def`` keyword.

    .. seealso:: :class:`Term` :class:`Typedef` :class:`Instance`

    """
    def __init__(self, stanza_type="Term", **kwargs):
        """
        Initialize from keyword arguments.
        """
        self.stanza_type = stanza_type
        self.tag_values = []
        self.values = {}

        sorted_tags = sorted(
            kwargs.iteritems(),
            key=lambda key_val: chr(1) if key_val[0] == "id" else key_val[0]
        )

        for tag, value in sorted_tags:
            if isinstance(value, basestring):
                tag, value, modifiers, comment = \
                    parse_tag_value(name_demangle(tag) + ": " + value)
            elif isinstance(value, tuple):
                tag, value, modifiers, comment = \
                    ((name_demangle(tag),) + value + (None, None))[:4]
            self.add_tag(tag, value, modifiers, comment)

        self.related = set()

    @property
    def is_annonymous(self):
        """
        Is this object anonymous.
        """
        value = self.get_values("is_annonymous")
        return bool(value)

    @property
    def id(self):
        """
        The id of this object.
        """
        value = self.get_values("id")
        return value[0] if value else None

    @property
    def name(self):
        """
        Name of this object
        """
        value = self.get_values("name")
        return value[0] if value else None

    def name_mangle(self, tag):
        return name_mangle(tag)

    def name_demangle(self, tag):
        return name_demangle(tag)

    def add_tag(self, tag, value, modifiers=None, comment=None):
        """
        Add `tag`, `value` pair to the object with optional modifiers and
        comment.

        Example::

            >>> term = OBOObject("Term")
            >>> term.add_tag("id", "FOO:002", comment="This is an id")
            >>> print term
            [Term]
            id: FOO:002 ! This is an id

        """
        tag = intern(tag)  # a small speed and memory benefit
        self.tag_values.append((tag, value, modifiers, comment))
        self.values.setdefault(tag, []).append(value)

    def add_tags(self, tag_value_iter):
        for tag, value, modifiers, comment in tag_value_iter:
            self.tag_values.append((tag, value, modifiers, comment))
            self.values.setdefault(tag, []).append(value)

    def update(self, other):
        """
        Update the term with tag value pairs from `other`
        (:class:`OBOObject`). The tag value pairs are appended to the
        end except for the `id` tag.

        """
        for tag, value, modifiers, comment in other.tag_values:
            if tag != "id":
                self.add_tag(tag, value, modifiers, comment)

    def get_values(self, tag):
        try:
            return self.values[tag]
        except KeyError:
            return []

    def tag_count(self):
        """
        Return the number of tags in this object.
        """
        return len(self.tag_values)

    def tags(self):
        """
        Return an list of all (tag, value, modifiers, comment) tuples.
        """
        return list(self.tag_values)

    def _format_single_tag(self, index):
        """
        Return a formated string representing index-th tag pair value.

        Example::

            >>> term = OBOObject(
            ...     "Term", id="FOO:001", name="bar",
            ...      def_="Example definition {modifier=frob} ! Comment")
            ...
            >>> term._format_single_tag(0)
            'id: FOO:001'
            >>> term._format_single_tag(1)
            'def: Example definition { modifier=frob } ! Comment'

        ..
            Why by index, and not by tag? Multiple tags are allowed.

        """
        tag, value, modifiers, comment = self.tag_values[index]
        res = ["%s: %s" % (tag, value)]
        if modifiers:
            res.append("{ %s }" % modifiers)
        if comment:
            res.append("! " + comment)
        return " ".join(res)

    def format_stanza(self):
        """
        Return a string stanza representation of this object.
        """
        stanza = ["[%s]" % self.stanza_type]
        for i in range(self.tag_count()):
            stanza.append(self._format_single_tag(i))
        return "\n".join(stanza)

    @classmethod
    def parse_stanza(cls, stanza):
        r'''
        Parse and return an OBOObject instance from a stanza string.

        >>> term = OBOObject.parse_stanza("""\
        ... [Term]
        ... id: FOO:001
        ... name: bar
        ... """)
        >>> print term.id, term.name
        FOO:001 bar

        '''
        lines = stanza.splitlines()
        stanza_type = lines[0].strip("[]")

        tag_values = [parse_tag_value(line) for line in lines[1:]
                      if ":" in line]

        obo = OBOObject(stanza_type)
        obo.add_tags(tag_values)
        return obo

    def related_objects(self):
        """
        Return a list of tuple pairs where the first element is
        relationship (typedef id) and the second object id whom the
        relationship applies to.

        """
        result = [(type_id, id)
                  for type_id in ["is_a"]  # TODO add other defined Typedef ids
                  for id in self.values.get(type_id, [])]

        result = result + [tuple(r.split(None, 1))
                           for r in self.values.get("relationship", [])]
        return result

    def __str__(self):
        """
        Return a string representation of the object in OBO format
        """
        return self.format_stanza()

    def __repr__(self):
        return ("{0.__name__}(id={1.id!r}, name={1.name}, ...)"
                .format(type(self), self))

    def __iter__(self):
        """
        Iterates over sub terms
        """
        return iter(self.related_objects())


class Term(OBOObject):
    """
    A 'Term' object in the ontology.
    """
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Term", *args, **kwargs)


class Typedef(OBOObject):
    """
    A 'Typedef' object in the ontology.
    """
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Typedef", *args, **kwargs)


class Instance(OBOObject):
    """
    An 'Instance' object in the ontology
    """
    def __init__(self, *args, **kwargs):
        OBOObject.__init__(self, "Instance", *args, **kwargs)


class OBOParser(object):
    r''' A simple parser for .obo files (inspired by xml.dom.pulldom)

    >>> from StringIO import StringIO
    >>> file = StringIO("""\
    ... header_tag: header_value
    ... [Term]
    ... id: FOO:001 { modifier=bar } ! comment
    ... """)
    >>> parser = OBOParser(file)
    >>> for event, value in parser:
    ...     print event, value
    ...
    HEADER_TAG ['header_tag', 'header_value']
    START_STANZA Term
    TAG_VALUE ('id', 'FOO:001', 'modifier=bar', 'comment')
    CLOSE_STANZA None

    '''
    def __init__(self, file):
        self.file = file

    def parse(self, progress_callback=None):
        """
        Parse the file and yield parse events.

        .. todo List events and values

        """
        data = self.file.read()
        header = data[: data.index("\n[")]
        body = data[data.index("\n[") + 1:]
        for line in header.splitlines():
            if line.strip():
                yield "HEADER_TAG", line.split(": ", 1)

        current = None
        #  For speed make these functions local
        startswith = str.startswith
        endswith = str.endswith
        parse_tag_value_ = parse_tag_value

        for line in body.splitlines():
            if startswith(line, "[") and endswith(line, "]"):
                yield "START_STANZA", line.strip("[]")
                current = line
            elif startswith(line, "!"):
                yield "COMMENT", line[1:]
            elif line:
                yield "TAG_VALUE", parse_tag_value_(line)
            else:  # empty line is the end of a term
                yield "CLOSE_STANZA", None
                current = None
        if current is not None:
            yield "CLOSE_STANZA", None

    def __iter__(self):
        """
        Iterate over parse events (same as parse())
        """
        return self.parse()


class OBOOntology(object):
    """
    An class representing an OBO ontology.

    :param file-like file:
        A optional file like object describing the ontology in obo format.

    """

    BUILTINS = BUILTIN_OBO_OBJECTS

    def __init__(self, file=None):
        self.objects = []
        self.header_tags = []
        self.id2term = {}
        self.alt2id = {}
        self._resolved_imports = []
        self._invalid_cache_flag = False
        self._related_to = {}

        # First load the built in OBO objects
        builtins = StringIO("\n" + "\n\n".join(self.BUILTINS) + "\n")
        self.load(builtins)
        if file:
            self.load(file)

    def add_object(self, obj):
        """
        Add an :class:`OBOObject` instance to this ontology.
        """
        if obj.id in self.id2term:
            raise ValueError("OBOObject with id: %s already in "
                             "the ontology" % obj.id)
        self.objects.append(obj)
        self.id2term[obj.id] = obj
        self._invalid_cache_flag = True

    def add_header_tag(self, tag, value):
        """
        Add header tag, value pair to this ontology.
        """
        self.header_tags.append((tag, value))

    def load(self, file, progress_callback=None):
        """
        Load terms from a file.

        :param file-like file:
            A file-like like object describing the ontology in obo format.
        :param function progress_callback:
            An optional function callback to report on the progress.

        """
        if isinstance(file, basestring):
            file = open(file, "rb")

        parser = OBOParser(file)
        current = None
        tag_values = []
        for event, value in parser.parse(progress_callback=progress_callback):
            if event == "TAG_VALUE":
                tag_values.append(value)
            elif event == "START_STANZA":
                current = OBOObject(value)
            elif event == "CLOSE_STANZA":
                current.add_tags(tag_values)
                self.add_object(current)
                current = None
                tag_values = []
            elif event == "HEADER_TAG":
                self.add_header_tag(*value)
            elif event != "COMMENT":
                raise Exception("Parse Error! Unknown parse "
                                "event {0}".format(event))

        imports = [value for tag, value in self.header_tags
                   if tag == "import"]

        if imports:
            warnings.warn("Import header tags are not supported")

#        while imports:
#            url = imports.pop(0)
#            if uri not in self._resolved_imports:
#                imported = self.parse_file(open(url, "rb"))
#                ontology.update(imported)
#                self._resolved_imports.append(uri)

    def dump(self, file):
        # deprecated use write
        self.write(file)

    def write(self, stream):
        """
        Write the contents of the ontology to a `file` in .obo format.

        :param file-like file:
            A file like object.

        """
        if isinstance(stream, basestring):
            stream = open(stream, "wb")

        for key, value in self.header_tags:
            stream.write(key + ": " + value + "\n")

        # Skip the builtins
        for obj in self.objects[len(self.BUILTINS):]:
            stream.write("\n")
            stream.write(obj.format_stanza())
            stream.write("\n")

    def update(self, other):
        """
        Update this ontology with the terms from `other`.
        """
        for term in other:
            if term.id in self:
                if not term.is_annonymous:
                    self.term(term.id).update(term)
                else:  # Do nothing
                    pass
            else:
                self.add_object(term)
        self._invalid_cache_flag = True

    def _cache_validate(self, force=False):
        """
        Update the relations cache if `self._invalid_cache` flag is set.
        """
        if self._invalid_cache_flag or force:
            self._cache_relations()

    def _cache_relations(self):
        """
        Collect all relations from parent to a child and store it in
        ``self._related_to`` member.

        """
        related_to = defaultdict(list)
        for obj in self.objects:
            for rel_type, id in self.related_terms(obj):
                term = self.term(id)
                related_to[term].append((rel_type, obj))

        self._related_to = related_to
        self._invalid_cache_flag = False

    def term(self, id):
        """
        Return the :class:`OBOObject` associated with this id.

        :param str id:
            Term id string.

        """
        if isinstance(id, basestring):
            if id in self.id2term:
                return self.id2term[id]
            elif id in self.alt2id:
                return self.id2term[self.alt2id[id]]
            else:
                raise ValueError("Unknown term id: %r" % id)
        elif isinstance(id, OBOObject):
            return id

    def terms(self):
        """
        Return all :class:`Term` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanza_type == "Term"]

    def term_by_name(self, name):
        """
        Return the term with name `name`.
        """
        terms = [t for t in self.terms() if t.name == name]
        if len(terms) != 1:
            raise ValueError("Unknown term name: %r" % name)
        return terms[0]

    def typedefs(self):
        """
        Return all :class:`Typedef` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanza_type == "Typedef"]

    def instances(self):
        """
        Return all :class:`Instance` instances in the ontology.
        """
        return [obj for obj in self.objects if obj.stanza_type == "Instance"]

    def root_terms(self):
        """
        Return all root terms (terms without any parents).
        """
        return [term for term in self.terms() if not self.parent_terms(term)]

    def related_terms(self, term):
        """
        Return a list of (`rel_type`, `term_id`) tuples where `rel_type` is
        relationship type (e.g. 'is_a', 'has_part', ...) and `term_id` is
        the id of the term in the relationship.

        """
        term = self.term(term) if not isinstance(term, OBOObject) else term
        related = [(tag, value)
                   for tag in ["is_a"]  # TODO: add other typedef ids
                   for value in term.values.get(tag, [])]
        relationships = term.values.get("relationship", [])
        for rel in relationships:
            related.append(tuple(rel.split(None, 1)))
        return related

    def edge_types(self):
        """
        Return a list of all edge types in the ontology.
        """
        return [obj.id for obj in self.objects if obj.stanza_type == "Typedef"]

    def parent_edges(self, term):
        """
        Return a list of (rel_type, parent_term) tuples.
        """
        term = self.term(term)
        parents = []
        for rel_type, parent in self.related_terms(term):
            parents.append((rel_type, self.term(parent)))
        return parents

    def child_edges(self, term):
        """
        Return a list of (rel_type, source_term) tuples.
        """
        self._cache_validate()
        term = self.term(term)
        return self._related_to.get(term, [])

    def super_terms(self, term):
        """
        Return a set of all super terms of `term` up to the most general one.
        """
        terms = self.parent_terms(term)
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(self.parent_terms(term) - visited)
        return visited

    def sub_terms(self, term):
        """
        Return a set of all sub terms for `term`.
        """
        terms = self.child_terms(term)
        visited = set()
        queue = set(terms)
        while queue:
            term = queue.pop()
            visited.add(term)
            queue.update(self.child_terms(term) - visited)
        return visited

    def child_terms(self, term):
        """
        Return a set of all child terms for this `term`.
        """
        self._cache_validate()
        term = self.term(term)
        children = []
        for rel_type, term in self.child_edges(term):
            children.append(term)
        return set(children)

    def parent_terms(self, term):
        """
        Return a set of all parent terms for this `term`.
        """
        term = self.term(term)
        parents = []
        for rel_type, id in self.parent_edges(term):
            parents.append(self.term(id))
        return set(parents)

    def relations(self):
        """
        Return a list of all relations in the ontology.
        """
        relations = []
        for obj in self.objects:
            for type_id, id in  obj.related:
                target_term = self.term(id)
            relations.append((obj, type_id, target_term))
        return relations

    def __len__(self):
        """
        Return the number of all objects in the ontology.
        """
        return len(self.objects)

    def __iter__(self):
        """
        Return an iterator over all objects in the ontology.
        """
        return iter(self.objects)

    def __contains__(self, oboid):
        return oboid in self.id2term

    def __getitem__(self, oboid):
        """
        Get the object by it's id `oboid`
        """
        return self.id2term[oboid]

    def to_network(self, terms=None):
        """
        Return an Orange.network.Network instance constructed from
        this ontology.

        """
        edge_types = self.edge_types()
        terms = self.terms()
        from Orange.orng import orngNetwork
        import orange

        network = orngNetwork.Network(len(terms), True, len(edge_types))
        network.objects = dict([(term.id, i) for i, term in enumerate(terms)])

        edges = defaultdict(set)
        for term in self.terms():
            related = self.related_terms(term)
            for relType, relTerm in related:
                edges[(term.id, relTerm)].add(relType)

        edgeitems = edges.items()
        for (src, dst), eTypes in edgeitems:
            network[src, dst] = [1 if e in eTypes else 0 for e in edge_types]

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
                                orange.EnumVariable("relationship", values=list(edge_types))
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

    def to_networkx(self, terms=None):
        """
        Return a NetworkX graph of this ontology.
        """
        import networkx
        graph = networkx.Graph()

        edge_types = self.edge_types()

        edge_colors = {"is_a": "red"}

        if terms is None:
            terms = self.terms()
        else:
            terms = [self.term(term) for term in terms]
            super_terms = [self.super_terms(term) for term in terms]
            terms = reduce(set.union, super_terms, set(terms))

        for term in terms:
            graph.add_node(term.id, name=term.name)

        for term in terms:
            for rel_type, rel_term in self.related_terms(term):
                rel_term = self.term(rel_term)
                if rel_term in terms:
                    graph.add_edge(term.id, rel_term.id, label=rel_type,
                                   color=edge_colors.get(rel_type, "blue"))

        return graph

    def to_graphviz(self, terms=None):
        """
        Return an pygraphviz.AGraph representation of the ontology.
        If `terms` is not `None` it must be a list of terms in the ontology.
        The graph will in this case contain only the super graph of those
        terms.

        """
        import pygraphviz as pgv
        graph = pgv.AGraph(directed=True, name="ontology")

        edge_types = self.edge_types()

        edge_colors = {"is_a": "red"}

        if terms is None:
            terms = self.terms()
        else:
            terms = [self.term(term) for term in terms]
            super_terms = [self.super_terms(term) for term in terms]
            terms = reduce(set.union, super_terms, set(terms))

        for term in terms:
            graph.add_node(term.id, label=term.name)

        for root in self.root_terms():
            node = graph.get_node(root.id)
            node.attr["rank"] = "max"

        for term in terms:
            for rel_type, rel_term in self.related_terms(term):
                rel_term = self.term(rel_term)
                if rel_term in terms:
                    graph.add_edge(term.id, rel_term.id, label=rel_type,
                                   color=edge_colors.get(rel_type, "blue"))

        return graph


def name_mangle(tag):
    """
    Mangle tag name if it conflicts with python keyword.

    >>> term.name_mangle("def"), term.name_mangle("class")
    ('def_', 'class_')

    """
    if keyword.iskeyword(tag):
        return tag + "_"
    else:
        return tag


def name_demangle(tag):
    """
    Reverse of `name_mangle`.
    """
    if tag.endswith("_") and keyword.iskeyword(tag[:-1]):
        return tag[:-1]
    else:
        return tag


def load(file):
    """
    Load an ontology from a .obo file.
    """
    return OBOOntology(file)


def foundry_ontologies():
    """
    Return a list of ontologies available from the OBOFoundry `website
    <http://www.obofoundry.org/>`_. The list contains a tuples of the form
    `(title, url)` for instance
    ``('Biological process', 'http://purl.obolibrary.org/obo/go.obo')``

    """
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
    seinfeld = StringIO("""
[Typedef]
id: parent

[Typedef]
id: child
inverse_of: parent ! not actually used yet

[Term]
id: 001
name: George

[Term]
id: 002
name: Estelle
relationship: parent 001 ! George

[Term]
id: 003
name: Frank
relationship: parent 001 ! George

""")  # TODO: fill the ontology with all characters
    term = OBOObject.parse_stanza(stanza)

    seinfeld = OBOOntology(seinfeld)
    print seinfeld.child_edges("001")

    doctest.testmod(extraglobs={"stanza": stanza, "term": term},
                    optionflags=doctest.ELLIPSIS)
