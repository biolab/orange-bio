import orange
#import mMisc as m
import numpy
import stats
#import mOrngData
import random

"""
Gene set enrichment analysis.

Author: Marko Toplak
"""

"""
Correlation methods.
"""

def mean(l):
    return float(sum(l))/len(l)

def rankingFromOrangeMeas(meas):
    """
    Creates a function that sequentally ranks all attributes and returns
    results in a list. Ranking function is build out of 
    orange.MeasureAttribute.
    """
    return lambda d: [ meas(i,d) for i in range(len(data.domain.attributes)) ]

class MA_pearsonCorrelation:
    """
    Calling an object of this class computes Pearson correlation of all
    attributes against class.
    """

    def __call__(self, i, data):
        dom2 = orange.Domain([data.domain.attributes[i]], data.domain.classVar)
        data2 = orange.ExampleTable(dom2, data)
        a,c = data2.toNumpy("A/C")
        return numpy.corrcoef(c,a[:,0])[0,1]

class MA_signalToNoise:
    """
    Returns signal to noise measurement: difference of means of two classes
    divided by the sum of standard deviations for both classes. 
    """

    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b

    def __call__(self, i, data):
        cv = data.domain.classVar

        #for faster computation. to save dragging many attributes along
        dom2 = orange.Domain([data.domain.attributes[i]], data.domain.classVar)
        data = orange.ExampleTable(dom2, data)
        i = 0

        if self.a == None: self.a = cv.values[0]
        if self.b == None: self.b = cv.values[1]

        def stdev(l):
            return stats.stdev(l)
    
        def stdevm(l):
            m = mean(l)
            std = stdev(l)
            #print std, 0.2*abs(1.0 if m == 0 else m)
            #return minmally 2*|mi|, where mi=0 is adjusted to mi=1
            return max(std, 0.2*abs(1.0 if m == 0 else m))

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[cv] == value]
    
        exa = avWCVal(self.a)
        exb = avWCVal(self.b)

        rval = (mean(exa)-mean(exb))/(stdevm(exa)+stdevm(exb))

        return rval



def enrichmentScoreRanked(subset, lcor, p=1.0):
    """
    Input data and subset. 
    
    subset: list of attribute indices of the input data belonging
        to the same set.
    lcor: correlations with class for each attribute in a list. 

    Returns enrichment score on given data.
    """

    ordered = [ (i,a) for i,a in enumerate(lcor) ] #original pos + correlation
    ordered.sort(lambda x,y: cmp(y[1],x[1])) #sort by correlation, descending
    ordered = nth(ordered, 0) #contains positions in the original list

    #add if gene is not in the subset
    notInA = -(1. / (len(lcor)-len(subset)))
    #base for addition if gene is in the subset
    inAb = 1./sum([abs(lcor[i])**p for i in subset])

    ess = [0.0]
    for i in ordered:
        ess.append(ess[-1] + \
            (inAb*abs(lcor[i])**p if i in subset else notInA)
        )

    maxEs = max(ess)
    minEs = min(ess)
    
    return (maxEs if abs(maxEs) > abs(minEs) else minEs, ess[1:])

#from mOrngData
def shuffleAttribute(data, attribute, locations):
    """
    Destructive!
    """
    attribute = data.domain[attribute]
    l = [None]*len(data)
    for i in range(len(data)):
        l[locations[i]] = data[i][attribute]
    for i in range(len(data)):
        data[i][attribute] = l[i]

def shuffleClass(data, rand=random.Random(0)):
    """
    Returns a dataset with values of class attribute randomly shuffled.
    """
    d2 = orange.ExampleTable(data.domain, data)
    locations = range(len(data))
    rand.shuffle(locations)
    shuffleAttribute(d2, d2.domain.classVar, locations)
    return d2

def shuffleList(l, rand=random.Random(0)):
    """
    Returns a copy of a shuffled input list.
    """
    import copy
    l2 = copy.copy(l)
    rand.shuffle(l2)
    return l2

def shuffleAttributes(data, rand=random.Random(0)):
    """
    Returns a dataset with a new attribute order.
    """
    natts = shuffleList(list(data.domain.attributes), rand)
    dom2 = orange.Domain(natts, data.domain.classVar)
    d2 = orange.ExampleTable(dom2, data)
    return d2

def gseapval(es, esnull):
    """
    From article (PNAS):
    estimate nominal p-value for S from esnull by using the positive
    or negative portion of the distribution corresponding to the sign 
    of the observed ES(S).

    WHAT DOES THAT MEAN?
    """
    
    if es < 0:
        return float(len([ a for a in esnull if a <= es ]))/ \
            len([ a for a in esnull if a < 0])    
    else: 
        return float(len([ a for a in esnull if a >= es ]))/ \
            len([ a for a in esnull if a > 0])


def enrichmentScore(data, subset, rankingf):
    """
    Returns enrichment score and running enrichment score.
    """
    lcor = rankingf(data)
    es,l = enrichmentScoreRanked(subset, lcor)
    return es,l

def gseaE(data, subsets, rankingf=rankingFromOrangeMeas(MA_signalToNoise()), \
        n=100, permutation="class", *kwargs):
    """
    Run GSEA algorithm on an example table.

    data: orange example table. 
    subsets: list of distinct subsets of data.
    rankingf: function that returns correlation to class of each 
        variable.
    n: number of random permutations to sample null distribution.
    permutation: "class" for permutating class, else permutate attribute 
        order.

    """
    enrichmentScores = []
 
    for subset in subsets:

        es = enrichmentScore(data, subset, rankingf)[0]
        enrichmentScores.append(es)

    runOptCallbacks(kwargs)

    #save correlations to speed up attribute order permutation
    lcor = []
    if permutation != "class":
        lcor = rankingf(data)

    enrichmentNulls = [ [] for a in range(len(subsets)) ]

    for i in range(n):

        if permutation == "class":
            d2 = shuffleClass(data, random.Random(2000+i)) #fixed permutation
        else:
            #SPEEDUP by saving correlations and skipping their recalculation
            r2 = shuffleList(lcor, random.Random(2000+i))

        for si,subset in enumerate(subsets):

            #type of permutation
            if permutation == "class":
                esn = enrichmentScore(d2, subset, rankingf)[0]
            else:
                esn = enrichmentScoreRanked(subset, r2)[0]

            enrichmentNulls[si].append(esn)

        runOptCallbacks(kwargs)


    return gseaSignificance(enrichmentScores, enrichmentNulls)


def runOptCallbacks(rargs):
    if "callback" in rargs:
        try:
            [ a() for a in rargs["callback"] ]
        except:
            rargs["callback"]()
            
        

def gseaR(rankings, subsets, n=100, **kwargs):
    """
    """

    enrichmentScores = []
 
    for subset in subsets:

        es = enrichmentScoreRanked(subset, rankings)[0]
        enrichmentScores.append(es)
    
    runOptCallbacks(kwargs)

    enrichmentNulls = [ [] for a in range(len(subsets)) ]

    for i in range(n):
        
        r2 = shuffleList(rankings, random.Random(2000+i))

        for si,subset in enumerate(subsets):

            esn = enrichmentScoreRanked(subset, r2)[0]
            enrichmentNulls[si].append(esn)

        runOptCallbacks(kwargs)

    return gseaSignificance(enrichmentScores, enrichmentNulls)


def gseaSignificance(enrichmentScores, enrichmentNulls):

    enrichmentPVals = []
    nEnrichmentScores = []
    nEnrichmentNulls = []

    for i in range(len(enrichmentScores)):
        es = enrichmentScores[i]
        enrNull = enrichmentNulls[i]

        enrichmentPVals.append(gseapval(es, enrNull))


        #normalize the ES(S,pi) and the observed ES(S), separetely rescaling
        #the positive and negative scores by divident by the mean of the 
        #ES(S,pi)

        meanPos = mean([a for a in enrNull if a >= 0])
        meanNeg = mean([a for a in enrNull if a < 0])

        def normalize(s):
            if s >= 0:
                return s/meanPos
            else:
                return -s/meanNeg

        nes = normalize(es)
        nEnrichmentScores.append(nes)
        
        nenrNull = [ normalize(s) for s in enrNull ]
        nEnrichmentNulls.append(nenrNull)
 
    #create a histogram of all NES(S,pi) over all S and pi
    vals = reduce(lambda x,y: x+y, nEnrichmentNulls)
    
    """
    Use this null distribution to compute an FDR q value, for a given NES(S) =
    NES* >= 0. The FDR is the ratio of the percantage of all (S,pi) with
    NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
    observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S)
    = NES* <= 0.
    """

    fdrs = []
    pvalsF = []

    import operator

    for i in range(len(enrichmentScores)):

        nes = nEnrichmentScores[i]
        if nes >= 0:
            op0 = operator.ge
            opn = operator.ge
        else:
            op0 = operator.lt
            opn = operator.le

        allPos = [a for a in vals if op0(a,0)]
        allHigherAndPos = [a for a in allPos if opn(a,nes) ]
        top = len(allHigherAndPos)/float(len(allPos)) #p value

        pvalsF.append(top)

        nesPos = [a for a in nEnrichmentScores if op0(a,0) ]
        nesHigherAndPos = [a for a in nesPos if opn(a,nes) ]
        down = len(nesHigherAndPos)/float(len(nesPos))

        fdrs.append(top/down)
    
    return zip(enrichmentScores, nEnrichmentScores, enrichmentPVals, fdrs)

import orngKEGG
#orngKEGG.default_database_path = "/home/marko/cached/data/kegg/"

def nth(l,n): return [ a[n] for a in l ]

def inverseDic(dic):
    item = dic.items()
    return dict(zip(nth(item,1),nth(item,0)))

class GeneMatcher(object):

    def __init__(self, caseSensitive=False):
        pass

    def setTargets(self, targets):
        pass

    def match(self, genes):
        pass

class GSEA(object):

    def __init__(self, organism="hsa"):
        self.a = "blabla"
        self.trans = {}
        self.keggorg = orngKEGG.KEGGOrganism(organism) #to do translation
        self.genesets = {}
        self.attrnames = []

    def addTransl(self, trans):
        self.trans.update(trans)
        
    def _translate(self, a):
        if a not in self.trans:
            uid,_,_ = self.keggorg.get_unique_gene_ids([a], caseSensitive=True) #unique id
            if len(uid) > 0:
                return uid[uid.keys()[0]]
            else:
                return None
        else:
            return self.trans[a]

    def translate(self, a):
        return (a.lower(), self._translate(a.lower())) #both a.lower !!

    def matchTargets(self):
        td = [ self.translate(a) for a in self.attrnames]

        def leaveOne(x):
            if x[1] == None:
                return x[0]
            else:
                return x[1]

        return map(leaveOne, td)

    def keepOnlyMeanAttrs(self, data):

        def attrOk(a):
            if a.varType != orange.VarTypes.Continuous:
                return False
            dom2 = orange.Domain([a],False)
            d2 = orange.ExampleTable(dom2, data)
            vals = [ex[0].value for ex in d2 if not ex[0].isSpecial()]
            if len(vals) < 1:
                return False
            #TODO. which attributes can I cancel with larger tables
            return True
            
        natts = [ a for a in data.domain.attributes if attrOk(a) ]
        ndom = orange.Domain(natts, data.domain.classVar)

        return orange.ExampleTable(ndom, data)

    def setData(self, data):

        data = self.keepOnlyMeanAttrs(data)

        self.data = data

        attrnames = [ a.name for a in data.domain.attributes ]

        #to from lowercase names to real names
        trn = {}
        for a in attrnames:
            trn[a.lower()] = a #first a.lower !!
        self.toRealNames = trn
        
        translations = [ self.translate(a) for a in attrnames ]
        translations = filter(lambda x: x[1] != None, translations)

        self.addTransl(translations)

        self.attrnames = attrnames
        self.targets = self.matchTargets()

    def genecompare(self, gene):
        gene = gene.lower() #loweeq
        transl =  self.translate(gene)[1]
        if transl != None:
            return transl
        else:
            return gene

    def dataMatch(self, genes):
        """
        Function returns a dictionary of an old value: matching data
        """
        targets = self.targets
        #print targets

        targetmap = dict(zip(targets,[1]*len(targets)))
       
        def matchingTarget(gene):
            """
            Find a match in input data for a given gene.
            """
            if self.genecompare(gene) in targetmap:
                return self.genecompare(gene)
            else:
                return None

        matches = [ (gene,matchingTarget(gene)) for gene in genes if matchingTarget(gene) ]

        def reverse(gene):
            id = inverseDic(self.trans)
            if gene in id:
                return id[gene]
            else:
                return gene

        matches = [ (a,reverse(b)) for a,b in matches ]

        return matches

    class AddGenesetException(Exception): pass

    def addGeneset(self, genesetname, genes):
        if genesetname in self.genesets:
            raise GenesetNameException("Geneset with the name " + \
                + genesetname + " is already in genesets.")
        else:
            #print genesetname, genes
            #matching to unified gene names
            datamatch = self.dataMatch(genes)
    
            #print "Matched", len(datamatch), "of", len(genes)

            #calulate coverage
            self.genesets[genesetname] = ( genes, datamatch )

            #print self.genesets

    def compute(self, minSize=3, maxSize=1000, minPart=0.1, n=100, **kwargs):

        subsets = self.genesets.items()

        data = self.data

        def okSizes(orig, transl):
            if len(transl) >= minSize and len(transl) <= maxSize \
                and float(len(transl))/len(orig) >= minPart:
                return True
            return False

        subsetsok = [ (a,(b,c)) for a,(b,c) in subsets if okSizes(b,c) ]
        
        #we need a mapping from attribute names to their

        namesToIndices = dict( \
            [(at.name, i) for i,at in enumerate(data.domain.attributes)])

        #print namesToIndices
        #print self.toRealNames

        nsubsets = []
        nsubsetsNames = []

        for subset in subsetsok:
            nsubsets.append( \
                [namesToIndices[self.toRealNames[b]] for a,b in subset[1][1]])
            nsubsetsNames.append([self.toRealNames[b] for a,b in subset[1][1]])

        #print nsubsets

        if len(self.data) > 1:
            gseal = gseaE(data, nsubsets, n=n, **kwargs)
        else:
            rankings = [ data[0][at].native() for at in data.domain.attributes ]
            gseal = gseaR(rankings, nsubsets, n=n, **kwargs)

        res = {}
    
        for i,subset in enumerate(subsetsok):
            name = subset[0]
            oSize = len(subset[1][0])
            tSize = len(subset[1][1])
            res[name] = list(gseal[i]) + [oSize, tSize] + [nsubsetsNames[i]]

        return res

if  __name__=="__main__":

    import mOrngData

    data = orange.ExampleTable("iris.tab")
    data = orange.ExampleTable("DLBCL.tab")
    data = mOrngData.getAttributes(data, range(100))

    subsets = [[2,3], [1,4,5,6], [3,4,5], [7,8,9], range(10,20), [14,14,20], [45,51,64], [33,46,49], [31,23,66], [3,54,55], [43,12,96], [1,34,12], [1,43,21], [5,34,87]]

    gseal = gseaE(data, subsets, rankingFromOrangeMeas(MA_signalToNoise()), n=10, permutation="class")

    #b =  rankingFromOrangeMeas(MA_signalToNoise())(data)
    #gseal = gseaR(b, subsets, n=1000)

    #for el in gseal:
    #    print el
 
    print gseal
