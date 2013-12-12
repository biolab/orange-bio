def only_option(a):
    if len(a) == 1:
        return list(a)[0]
    else:
        raise Exception()

class GenesetRegException(Exception): pass

class GeneSet(object):
    """ A single set of genes.
    """

    def __init__(self, genes=[], name=None, id=None, \
        description=None, link=None, organism=None, hierarchy=None, pair=None):
        """
        :param pair: Backward compatibility: convert a tuple (name, genes)
            into this object.
        """

        self.hierarchy = hierarchy     
        """ Hierarchy should be formated as a tuple, for example ``("GO", "biological_process")``"""

        self.genes = set(genes)
        """ A set of genes. Genes are strings. """

        self.name = name
        """ Gene set name. """

        self.id = id
        """ Short gene set ID. """

        self.description = description
        """ Gene set description. """

        self.link = link
        """ Link to further information about this gene set. """

        self.organism = organism
        """ Organism as a NCBI taxonomy ID. """

        if pair:
            self.id, self.genes = pair[0], set(pair[1])

    """
    the following functions are needed for sets of gene sets to be able
    to assess equality
    """

    def __hash__(self):
        return self.id.__hash__() + self.name.__hash__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def size(self):
        return len(self.genes)

    def cname(self, source=True, name=True):
        """ Return a gene set name with hieararchy. """
        oname = self.id
        if source and self.hierarchy:
            oname = "[ " + ", ".join(self.hierarchy) + " ] " + oname
        if name and self.name:
            oname = oname + " " + self.name
        return oname

    def to_odict(self, source=True, name=True):
        """
        Backward compatibility: returns a gene set as a tuple
        (id, list of genes).
        """
        return self.cname(source=source, name=name), self.genes

    def __repr__(self):
        return "GeneSet(" + ", ".join( [ 
            "id=" + str(self.id),
            "genes=" + str(self.genes),
            "name=" + str(self.name),
            "link=" + str(self.link),
            "hierarchy=" + str(self.hierarchy)
        ]) + ")"

class GeneSetIDException(Exception):
    pass

class GeneSets(set):
    """ A collection of gene sets: contains :class:`GeneSet` objects. 
    It is a subclass of Python's :obj:`set`. 
    """
    
    def __init__(self, input=None):
        """
        If input is a dictionary, the gene sets are converted to the current format.
        """
        if input != None and len(input) > 0:
            self.update(input)

    def update(self, input):
        if isinstance(input, GeneSets):
            super(GeneSets, self).update(input)
        else:
            prepared_genesets = [] #parse them all before adding,
                                   #so that it fails on error
            if hasattr(input, "items"):
                for i, g in input.items():
                    prepared_genesets.append(GeneSet(pair=(i, g)))
            else:
                for i in input:
                    if isinstance(i, GeneSet):
                        prepared_genesets.append(i)
                    else:
                        i, g = i
                        prepared_genesets.append(GeneSet(pair=(i, g)))

            for g in prepared_genesets:
                self.add(g)

    def to_odict(self):
        """ Return gene sets in old dictionary format. """
        return dict(gs.to_odict() for gs in self)

    def set_hierarchy(self, hierarchy):
        """ Sets hierarchy for all gene sets. """
        for gs in self:
            gs.hierarchy = hierarchy

    def __repr__(self):
        return "GeneSets(" + set.__repr__(self) + ")"

    def common_org(self):
        """ Return a common organism. """
        if len(self) == 0:
            raise GenesetRegException("Empty gene sets.")

        organisms = set(a.organism for a in self)

        try:
            return only_option(organisms)
        except:
            raise GenesetRegException("multiple organisms: " + str(organisms))

    def hierarchies(self):
        """ Return all hierarchies. """
        if len(self) == 0:
            raise GenesetRegException("Empty gene sets.")
        return set(a.hierarchy for a in self)

    def common_hierarchy(self):
        """ Return a common hierarchy. """
        hierarchies = self.hierarchies()

        def common_hierarchy1(hierarchies):
            def hier(l): return set(map(lambda x: x[:currentl], hierarchies))
            currentl = max(map(len, hierarchies))
            while len(hier(currentl)) > 1:
                currentl -= 1
            return only_option(hier(currentl))

        return common_hierarchy1(hierarchies)

    def split_by_hierarchy(self):
        """ Split gene sets by hierarchies. Return a list of :class:`GeneSets` objects. """
        hd = dict((h,GeneSets()) for h in  self.hierarchies())
        for gs in self:
            hd[gs.hierarchy].add(gs)
        return hd.values()

