import orngKEGG

def nth(l,n): return [ a[n] for a in l ]

def inverseDic(dic):
    item = dic.items()
    return dict(zip(nth(item,1),nth(item,0)))

class GeneMatcher(object):

    def __init__(self, targets, caseSensitive=True, organism="hsa"):

        self.trans = {}

        self.caseSensitive = caseSensitive

        self.keggorg = orngKEGG.KEGGOrganism(organism) #to do translation
        self.organism = organism

        attrnames = targets

        #from lowercase names to real names - if applicable. else identity
        trn = {}
        for a in attrnames:
            if self.caseSensitive:
                trn[a] = a
            else:
                trn[a.lower()] = a
        self.toRealNames = trn

        translations = [ self.translate(a) for a in attrnames ]
        translations = filter(lambda x: x[1] != None, translations)

        self.addTransl(translations)
        self.attrnames = attrnames

        self.targets = self.matchTargets()
        self.targetmap = dict(zip(self.targets,[1]*len(self.targets)))

    def addTransl(self, trans):
        self.trans.update(trans)
        self.transi = inverseDic(self.trans)

    def _translate(self, a):
        if a not in self.trans:
            uid,_,_ = self.keggorg.get_unique_gene_ids([a], caseSensitive=self.caseSensitive)
            if len(uid) > 0:
                return uid.keys()[0]
            else:
                return None
        else:
            return self.trans[a]

    def translate(self, a):
        if not self.caseSensitive:
            a = a.lower()
        return (a, self._translate(a))

    def matchTargets(self):
        td = [ self.translate(a) for a in self.attrnames]

        def leaveOne(x):
            if x[1] == None:
                return x[0]
            else:
                return x[1]

        return map(leaveOne, td)

    def genecompare(self, gene):
        if not self.caseSensitive:
            gene = gene.lower()
        transl =  self.translate(gene)[1]
        if transl != None:
            return transl
        else:
            return gene

    def match(self, genes):
        """
        Function returns a dictionary of an old value: matching data
        """
        targets = self.targets
        targetmap = self.targetmap

        def matchingTarget(gene):
            """
            Find a match in input data for a given gene.
            """
            gc = self.genecompare(gene)
            if gc in targetmap:
                return gc
            else:
                return None

        matches = [ (gene,matchingTarget(gene)) for gene in genes if matchingTarget(gene)]

        def reverse(gene):
            if gene in self.transi:
                return self.transi[gene]
            else:
                return gene

        matches = [ (a,self.toRealNames[reverse(b)]) for a,b in matches ]

        return matches

