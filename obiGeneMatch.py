import obiKEGG

def nth(l,n): return [ a[n] for a in l ]

def GeneMatch(targets, organism, caseSensitive=False):
    return MatcherSequence([MatchName(targets, caseSensitive=caseSensitive),
            MatchKEGG(targets, organism, caseSensitive=caseSensitive)])

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

def issequencens(x):
    "Is x a sequence and not string ? We say it is if it has a __getitem__ method and it is not an instance of basestring."
    return hasattr(x, '__getitem__') and not isinstance(x, basestring)

class Matcher(object):

    def targets(self, targets):
        notImplemented()

    def matchOne(self, gene):
        notImplemented()

    def matchL(self, genes):
        return [ (gene, self.matchOne(gene)) for gene in genes if self.matchOne(gene) ]

    def match(self, genes):
        """
        Genes can be a string or a single gene.
        """
        if issequencens(genes):
            return sorted(self.matchL(genes))
        else:
            return self.matchOne(genes)

class MatchName(Matcher):
    """
    Match genes by name. Mather matches genes with exactly the same
    name, possibly ignoring case. If case is not ignored, then this is really
    dumb "gene matching".
    """

    def __init__(self, targets, caseSensitive=False):
        self.caseSensitive = caseSensitive
        self.targets(targets)

    def targets(self, targets):
        self.toRealNames = constructToRealNames(targets, caseSensitive=self.caseSensitive)

    def matchOne(self, gene):
        if not self.caseSensitive:
            gene = gene.lower()
        return self.toRealNames.get(gene, None)

class MatchKEGG(Matcher):
    """
    Match by KEGG unique id.
    """
    def __init__(self, targets, organism, caseSensitive=False):
        self.caseSensitive = caseSensitive
        self.keggorg = obiKEGG.KEGGOrganism(organism)
        self.targets(targets)
    
    def targets(self, targets):
        uniqueids = [ (self.tounique(a),a) for a in targets ]
        uniqueids = filter(lambda x: x[0] != None, uniqueids)
        self.toname = dict(uniqueids)

    def tounique(self, a):
        if not self.caseSensitive:
            a = a.lower()

        uid,_,_ = self.keggorg.get_unique_gene_ids([a], caseSensitive=self.caseSensitive)

        if len(uid) > 0:
            return uid.keys()[0]
        else:
            return None

    def matchOne(self, gene):
        if not self.caseSensitive:
            gene = gene.lower()
        return self.toname.get(self.tounique(gene), None)

class MatcherSequence(Matcher):
    
    def __init__(self, matchers):
        self.matchers = matchers

    def matchOne(self, gene):
        for matcher in self.matchers:
            m = matcher.matchOne(gene)
            if m != None:
                return m
        return None

    def targets(self, targets):
        for matcher in self.matchers:
            matcher.targets(targets)


_dbOrgMap = {"goa_human":"hsa", "dictyBase":"ddi", "sgd":"sce", "fb":"dme", "tigr_Aphagocytophilum":"aph", "PAMGO_Atumefaciens":"atu", "tair":"ath", "tigr_Banthracis":"ban",
                "goa_cow":"bta", "tigr_Chydrogenoformans":"chy", "wb":"cel", "tigr_Cjejuni":"cje", "cgd":"cal", "tigr_Cperfringens":"cpf", "tigr_Cpsychrerythraea":"cps", "tigr_Cburnetii":"cbu",
                "zfin":"dre", "tigr_Dethenogenes":"det", "tigr_Echaffeensis":"ech", "goa_chicken":"gga", "tigr_Gsulfurreducens":"gsu", "tigr_Hneptunium":"hne", "GeneDB_Lmajor":"lma",
                "tigr_Lmonocytogenes":"lmf", "PAMGO_Mgrisea":"mgr", "tigr_Mcapsulatus":"mca", "mgi":"mmu", "tigr_Nsennetsu":"nse", "gramene_oryza":"osa", "GeneDB_Pfalciparum":"pfa",
                "pseudocap":"pae", "tigr_Pfluorescens":"pfl", "tigr_Psyringae":"pst", "tigr_Psyringae_phaseolicola":"psp", "rgd":"rno", "GeneDB_Spombe":"spo", "tigr_Soneidensis":"son",
                "tigr_Spomeroyi":"sil", "GeneDB_Tbrucei":"tbr", "tigr_Vcholerae":"vch", "ecocyc":"eco"}
class GeneMatchMk2(object):
    dbNameMap = {"UniProtKB":"UniProt", "SGD":"SGD", "dictyBase":"DictyBase"}
    dbOrgMap = _dbOrgMap
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
    attrnames =  [ str(a.name) for a in data.domain.attributes ] + [ "blabla1" ]
    print attrnames

    gm = GeneMatch(attrnames, organism="hsa", caseSensitive=False)
    #gm = MatchName(attrnames, caseSensitive=False)
    #gm = MatchKEGG(attrnames, "hsa", caseSensitive=False)

    #gm = MatcherSequence([MatchName([], caseSensitive=False), MatchKEGG([], "hsa", caseSensitive=False)])
    #gm.targets(attrnames)

    gen1 = getDefaultGenesets()

    gen1["TEST1"] = [ "Blabla1" ]

    t1 = time.time()
    res = []
    for k,gs in sorted(gen1.items()):
        res.append((k, gs, gm.match(gs)))
    t2 = time.time() - t1

    for a in res:
        print a[0], a[2]

    print "TIME", t2
