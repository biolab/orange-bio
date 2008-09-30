import obiKEGG

def nth(l,n): return [ a[n] for a in l ]

def inverseDic(dic):
    item = dic.items()
    return dict(zip(nth(item,1),nth(item,0)))

class GeneMatch(object):

    def __init__(self, targets, caseSensitive=True, organism="hsa"):

        self.trans = {}

        self.caseSensitive = caseSensitive

        self.keggorg = obiKEGG.KEGGOrganism(organism) #to do translation
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
        #print self.targets
        self.targetmap = dict(zip(self.targets,[1]*len(self.targets)))

    def addTransl(self, trans):
        self.trans.update(trans)
        self.transi = inverseDic(self.trans)

    def _translate(self, a):
        if a not in self.trans:
            uid,_,_ = self.keggorg.get_unique_gene_ids([a], caseSensitive=self.caseSensitive)
            if len(uid) > 0:
                #print a, uid
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
        #print td

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
        Function returns a tuple of an (old value, matching data)
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

        return sorted(matches)

def constructToRealNames(targets, caseSensitive):
    """
    Map to original gene name: from lowercase
    names to real names, if caseSensitive == False.
    Else it is an identity.
    """
    trn = {}
    for a in targets:
        if caseSensitive:
            trn[a] = a
        else:
            trn[a.lower()] = a
    return trn

class MatchName(object):
    """
    Match genes by name. Mather matches genes with exactly the same
    name, possibly ignoring case. If case is not ignored, then this is really
    dumb "gene matching". - not properly tested!
    """

    def __init__(self, targets, caseSensitive=False):

        self.caseSensitive = caseSensitive
        self.toRealNames = constructToRealNames(targets, caseSensitive=caseSensitive)
        print self.toRealNames

    def matchOne(self, gene):
        if not self.caseSensitive:
            gene = gene.lower()

        return self.toRealNames.get(gene, None)

    def matchL(self, genes):
        return [ (gene, self.matchOne(gene)) for gene in genes if self.matchOne(gene) ]

    def match(self, genes):
        return sorted(self.matchL(genes))

class MatchKEGG(object):
    """
    Match by KEGG unique id. - not tested properly!
    """
    def __init__(self, targets, caseSensitive=False, organism="hsa"):

        self.caseSensitive = caseSensitive

        self.keggorg = obiKEGG.KEGGOrganism(organism)

        attrnames = targets
        self.toRealNames = constructToRealNames(targets, caseSensitive=caseSensitive)

        uniqueids = [ (a, self.tounique(a)) for a in attrnames ]
        uniqueids = filter(lambda x: x[1] != None, translations)

        self.trans = {}

        self.addTransl(translations)
        self.attrnames = targets

        self.targets = self.matchTargets()
        #print self.targets
        self.targetmap = dict(zip(self.targets,[1]*len(self.targets)))

    def addTransl(self, trans):
        self.trans.update(trans)
        self.transi = inverseDic(self.trans)

    def tounique(self, a):
        if not self.caseSensitive:
            a = a.lower()

        uid,_,_ = self.keggorg.get_unique_gene_ids([a], caseSensitive=self.caseSensitive)

        if len(uid) > 0:
            return uid.keys()[0]
        else:
            return None

    def matchTargets(self):
        td = [ self.translate(a) for a in self.attrnames]
        #print td

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
        Function returns a tuple of an (old value, matching data)
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

        return sorted(matches)



class GeneMatchMk2(object):
    dbNameMap = {"UniProtKB":"UniProt", "SGD":"SGD", "dictyBase":"DictyBase"}
    dbOrgMap = {"goa_human":"hsa", "dictyBase":"ddi", "sgd":"sce", "fb":"dme", "tigr_Aphagocytophilum":"aph", "PAMGO_Atumefaciens":"atu", "tair":"ath", "tigr_Banthracis":"ban",
                "goa_cow":"cow", "tigr_Chydrogenoformans":"chy", "wb":"cel", "tigr_Cjejuni":"cje", "cgd":"cal", "tigr_Cperfringens":"cpf", "tigr_Cpsychrerythraea":"cps", "tigr_Cburnetii":"cbu",
                "zfin":"dre", "tigr_Dethenogenes":"det", "tigr_Echaffeensis":"ech", "goa_chicken":"gga", "tigr_Gsulfurreducens":"gsu", "tigr_Hneptunium":"hne", "GeneDB_Lmajor":"lma",
                "tigr_Lmonocytogenes":"lmf", "PAMGO_Mgrisea":"mgr", "tigr_Mcapsulatus":"mca", "mgi":"mmu", "tigr_Nsennetsu":"nse", "gramene_oryza":"osa", "GeneDB_Pfalciparum":"pfa",
                "pseudocap":"pae", "tigr_Pfluorescens":"pfl", "tigr_Psyringae":"pst", "tigr_Psyringae_phaseolicola":"psp", "rgd":"rno", "GeneDB_Spombe":"spo", "tigr_Soneidensis":"son",
                "tigr_Spomeroyi":"sil", "GeneDB_Tbrucei":"tbr", "tigr_Vcholerae":"vch"}
    def __init__(self, keggOrg, caseSensitive=True):
        self.keggOrg = keggOrg
        self.caseSensitive = caseSensitive

    def GetNamesFromDB(self, names, dbName):
        dbName = self.dbNameMap.get(dbName, dbName)
        k, c, u = self.keggOrg.get_unique_gene_ids(names, self.caseSensitive)
        mapper = {}
        for key, name in k.items():
            links = self.keggOrg.api._genes[self.keggOrg.org][key].get_db_links()
            if dbName in links:
                for link in links[db]:
                    mapper[link] = name
            else:
                u.append(name)
        return mapper, c, u

    def LinkWith(self, names, aliasMapper):
        mapper = dict([(name, aliasMapper[name]) for name in names if name in aliasMapper])
        names = [name for name in names if name not in aliasMapper]
        k, c, u = self.keggOrg.get_unique_gene_ids(names, self.caseSensitive)
        for key, name in k.items():
            altNames = self.keggOrg.api._genes[self.keggOrg.org][key].get_alt_names()
            nameSet = set([aliasMapper[altName] for altName in altNames if altName in aliasMapper])
            if len(nameSet)==1:
                mapper[name] = nameSet.pop()
            elif len(nameSet)==0:
                u.append(name)
            else:
                c.append(name)
        return mapped, c, u

if __name__ == "__main__":

    import time

    def getDefaultGenesets():

        def unpckGS(filename):
            import pickle
            f = open(filename,'rb')
            return pickle.load(f)

        import orngEnviron
        return unpckGS(orngEnviron.directoryNames["bufferDir"] + "/gsea/geneSets_MSIGDB.pck")

    import orange

    data = orange.ExampleTable("allData.tab")
    attrnames =  [ str(a.name) for a in data.domain.attributes ]
    print attrnames

    #gm = GeneMatch(attrnames, organism="hsa", caseSensitive=False)
    gm = MatchName(attrnames, caseSensitive=False)

    gen1 = getDefaultGenesets()

    t1 = time.time()
    res = []
    for k,gs in sorted(gen1.items()):
        res.append((k, gs, gm.match(gs)))

    t2 = time.time() - t1

    for a in res:
        print a[0], a[2]

    print "TIME", t2
