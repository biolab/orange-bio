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

        return matches

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

