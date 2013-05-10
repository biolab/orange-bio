def only_option(a):
    if len(a) == 1:
        return list(a)[0]
    else:
        raise Exception()

class GenesetRegException(Exception): pass

class GeneSet(object):

    def __init__(self, genes=None, name=None, id=None, \
        description=None, link=None, organism=None, hierarchy=None, pair=None):
        """
        pair can be (id, listofgenes) - it is used before anything else.
        """
        if genes == None:
            genes = []

        self.hierarchy = hierarchy     
        self.genes = set(genes)
        self.name = name
        self.id = id
        self.description = description
        self.link = link
        self.organism = organism

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
        """ Constructs a gene set name with the hierarchy. """
        oname = self.id
        if source and self.hierarchy:
            oname = "[ " + ", ".join(self.hierarchy) + " ] " + oname
        if name and self.name:
            oname = oname + " " + self.name
        return oname

    def to_odict(self, source=True, name=True):
        """
        Returns a pair (id, listofgenes), like in old format.
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
    
    def __init__(self, input=None):
        """
        odict are genesets in old dict format.
        gs are genesets in new format
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
        """ Sets hierarchy for all gene sets """
        for gs in self:
            gs.hierarchy = hierarchy

    def __repr__(self):
        return "GeneSets(" + set.__repr__(self) + ")"

    def common_org(self):
        """ Returns the common organism. """
        if len(self) == 0:
            raise GenesetRegException("Empty gene sets.")

        organisms = set(a.organism for a in self)

        try:
            return only_option(organisms)
        except:
            raise GenesetRegException("multiple organisms: " + str(organisms))

    def hierarchies(self):
        """ Returns all hierachies """
        if len(self) == 0:
            raise GenesetRegException("Empty gene sets.")
        return set(a.hierarchy for a in self)

    def common_hierarchy(self):
        hierarchies = self.hierarchies()

        def common_hierarchy1(hierarchies):
            def hier(l): return set(map(lambda x: x[:currentl], hierarchies))
            currentl = max(map(len, hierarchies))
            while len(hier(currentl)) > 1:
                currentl -= 1
            return only_option(hier(currentl))

        return common_hierarchy1(hierarchies)

    def split_by_hierarchy(self):
        """ Splits gene sets by hierarchies. """
        hd = dict((h,GeneSets()) for h in  self.hierarchies())
        for gs in self:
            hd[gs.hierarchy].add(gs)
        return hd.values()

