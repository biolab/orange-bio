import orange
#import mMisc as m
import numpy
import stats
#import mOrngData
import random
import time
import math, os
from obiExpression import *

"""
Gene set enrichment analysis.

Author: Marko Toplak
"""

collectionsPath = None
def iset(data):
    """
    Is data orange.ExampleTable?
    """
    return isinstance(data, orange.ExampleTable)

def issequencens(x):
    "Is x a sequence and not string ? We say it is if it has a __getitem__ method and is not string."
    return hasattr(x, '__getitem__') and not isinstance(x, basestring)

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
    return lambda d: [ meas(i,d) for i in range(len(d.domain.attributes)) ]

def orderedPointersCorr(lcor):
    """
    Return a list of integers: indexes in original
    lcor. Elements in the list are ordered by
    their lcor[i] value. Higher correlations first.
    """
    ordered = [ (i,a) for i,a in enumerate(lcor) ] #original pos + correlation
    ordered.sort(lambda x,y: cmp(y[1],x[1])) #sort by correlation, descending
    ordered = nth(ordered, 0) #contains positions in the original list
    return ordered

def enrichmentScoreRanked(subset, lcor, ordered, p=1.0, rev2=None):
    """
    Input data and subset. 
    
    subset: list of attribute indices of the input data belonging
        to the same set.
    lcor: correlations with class for each attribute in a list. 

    Returns enrichment score on given data.

    This implementation efficiently handles "sparse" genesets (that
    cover only a small subset of all genes in the dataset).
    """

    #print lcor

    subset = set(subset)

    if rev2 == None:
        def rev(l):
            return numpy.argsort(l)
        rev2 = rev(ordered)

    #add if gene is not in the subset
    notInA = -(1. / (len(lcor)-len(subset)))
    #base for addition if gene is in the subset
    cors = [ abs(lcor[i])**p for i in subset ]
    sumcors = sum(cors)

    #this should not happen
    if sumcors == 0.0:
        return (0.0, None)
    
    inAb = 1./sumcors

    ess = [0.0]
    
    map = {}
    for i in subset:
        orderedpos = rev2[i]
        map[orderedpos] = inAb*abs(lcor[i]**p)
        
    #print sorted(map.items())

    last = 0

    maxSum = minSum = csum = 0.0

    for a,b in sorted(map.items()):
        #print a, b
        #print csum
        diff = a-last
        csum += notInA*diff
        #print csum
        last = a+1
        
        if csum < minSum:
            minSum = csum
        
        csum += b

        if csum > maxSum:
            maxSum = csum

    #finish it
    diff = (len(ordered))-last
    csum += notInA*diff
    #print "LAST",csum
    if csum < minSum:
        minSum = csum

    #print "MY", (maxSum if abs(maxSum) > abs(minSum) else minSum)

    """
    #BY DEFINITION
    print "subset", subset

    for i in ordered:
        ess.append(ess[-1] + \
            (inAb*abs(lcor[i]**p) if i in subset else notInA)
        )
        if i in subset:
            print ess[-2], ess[-1]
            print i, (inAb*abs(lcor[i]**p))

    maxEs = max(ess)
    minEs = min(ess)
    
    print "REAL", (maxEs if abs(maxEs) > abs(minEs) else minEs, ess[1:])

    """
    return (maxSum if abs(maxSum) > abs(minSum) else minSum, [])

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

def shuffleClass(datai, rands=0):
    """
    Returns a dataset with values of class attribute randomly shuffled.
    If multiple dataset are on input shuffle them all with the same random seed.
    """
    def shuffleOne(data):
        rand = random.Random(rands)
        d2 = orange.ExampleTable(data.domain, data)
        locations = range(len(data))
        rand.shuffle(locations)
        shuffleAttribute(d2, d2.domain.classVar, locations)
        return d2

    if iset(datai):
        return shuffleOne(datai)
    else:
        return [ shuffleOne(data) for data in datai ]

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
    """
    
    try:
        if es < 0:
            return float(len([ a for a in esnull if a <= es ]))/ \
                len([ a for a in esnull if a < 0])    
        else: 
            return float(len([ a for a in esnull if a >= es ]))/ \
                len([ a for a in esnull if a >= 0])
    except:
        return 1.0


def enrichmentScore(data, subset, rankingf):
    """
    Returns enrichment score and running enrichment score.
    """
    lcor = rankingf(data)
    ordered = orderedPointersCorr(lcor)
    es,l = enrichmentScoreRanked(subset, lcor, ordered)
    return es,l

def gseaE(data, subsets, rankingf=None, \
        n=100, permutation="class", **kwargs):
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

    if not rankingf:
        rankingf=rankingFromOrangeMeas(MA_signalToNoise())

    enrichmentScores = []
 
    lcor = rankingf(data)
    #print lcor

    ordered = orderedPointersCorr(lcor)

    def rev(l):
        return numpy.argsort(l)

    rev2 = rev(ordered)

    for subset in subsets:
        es = enrichmentScoreRanked(subset, lcor, ordered, rev2=rev2)[0]
        enrichmentScores.append(es)

    runOptCallbacks(kwargs)

    #print "PERMUTATION", permutation

    enrichmentNulls = [ [] for a in range(len(subsets)) ]

    for i in range(n):

        if permutation == "class":
            d2 = shuffleClass(data, 2000+i) #fixed permutation
            r2 = rankingf(d2)

        else:
            r2 = shuffleList(lcor, random.Random(2000+i))

        ordered2 = orderedPointersCorr(r2)
        rev22 = rev(ordered2)
        for si,subset in enumerate(subsets):
            esn = enrichmentScoreRanked(subset, r2, ordered2, rev2=rev22)[0]
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

    if "permutation" in kwargs:
        if kwargs["permutation"] == "class":
            raise Exception("Only gene permutation possible")

    enrichmentScores = []
 
    ordered = orderedPointersCorr(rankings)
    
    def rev(l):
        return numpy.argsort(l)

    rev2 = rev(ordered)

    for subset in subsets:

        es = enrichmentScoreRanked(subset, rankings, ordered, rev2=rev2)[0]
        enrichmentScores.append(es)
    
    runOptCallbacks(kwargs)

    enrichmentNulls = [ [] for a in range(len(subsets)) ]

    for i in range(n):
        
        r2 = shuffleList(rankings, random.Random(2000+i))
        ordered2 = orderedPointersCorr(r2)
        rev22 = rev(ordered2)

        for si,subset in enumerate(subsets):

            esn = enrichmentScoreRanked(subset, r2, ordered2, rev2=rev22)[0]
            enrichmentNulls[si].append(esn)

        runOptCallbacks(kwargs)

    return gseaSignificance(enrichmentScores, enrichmentNulls)


def gseaSignificance(enrichmentScores, enrichmentNulls):

    #print enrichmentScores

    enrichmentPVals = []
    nEnrichmentScores = []
    nEnrichmentNulls = []

    for i in range(len(enrichmentScores)):
        es = enrichmentScores[i]
        enrNull = enrichmentNulls[i]
        #print es, enrNull

        enrichmentPVals.append(gseapval(es, enrNull))

        #normalize the ES(S,pi) and the observed ES(S), separetely rescaling
        #the positive and negative scores by divident by the mean of the 
        #ES(S,pi)

        #print es, enrNull

        def normalize(s):
            try:
                if s == 0:
                    return 0.0
                if s >= 0:
                    meanPos = mean([a for a in enrNull if a >= 0])
                    #print s, meanPos
                    return s/meanPos
                else:
                    meanNeg = mean([a for a in enrNull if a < 0])
                    #print s, meanNeg
                    return -s/meanNeg
            except:
                return 0.0 #return if according mean value is uncalculable


        nes = normalize(es)
        nEnrichmentScores.append(nes)
        
        nenrNull = [ normalize(s) for s in enrNull ]
        nEnrichmentNulls.append(nenrNull)
 

    #FDR computation
    #create a histogram of all NES(S,pi) over all S and pi
    vals = reduce(lambda x,y: x+y, nEnrichmentNulls, [])

    """
    Use this null distribution to compute an FDR q value, for a given NES(S) =
    NES* >= 0. The FDR is the ratio of the percantage of all (S,pi) with
    NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
    observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S)
    = NES* <= 0.
    """

    nvals = numpy.array(sorted(vals))
    nnes = numpy.array(sorted(nEnrichmentScores))

    fdrs = []

    import operator

    for i in range(len(enrichmentScores)):

        nes = nEnrichmentScores[i]

        """
        #Strighfoward but slow implementation follows in comments.
        #Useful as code description.
        
        
        if nes >= 0:
            op0 = operator.ge
            opn = operator.ge
        else:
            op0 = operator.lt
            opn = operator.le

        allPos = [a for a in vals if op0(a,0)]
        allHigherAndPos = [a for a in allPos if opn(a,nes) ]

        nesPos = [a for a in nEnrichmentScores if op0(a,0) ]
        nesHigherAndPos = [a for a in nesPos if opn(a,nes) ]

        top = len(allHigherAndPos)/float(len(allPos)) #p value
        down = len(nesHigherAndPos)/float(len(nesPos))
        
        l1 = [ len(allPos), len(allHigherAndPos), len(nesPos), len(nesHigherAndPos)]

        allPos = allHigherAndPos = nesPos =  nesHigherAndPos = 1

        """

        if nes >= 0:
            allPos = int(len(vals) - numpy.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(len(vals) - numpy.searchsorted(nvals, nes, side="left"))
            nesPos = len(nnes) - int(numpy.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = len(nnes) - int(numpy.searchsorted(nnes, nes, side="left"))
        else:
            allPos = int(numpy.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(numpy.searchsorted(nvals, nes, side="right"))
            nesPos = int(numpy.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = int(numpy.searchsorted(nnes, nes, side="right"))
           
        """
        #Comparing results
        l2 = [ allPos, allHigherAndPos, nesPos, nesHigherAndPos ]
        diffs = [ l1[i]-l2[i] for i in range(len(l1)) ]
        sumd = sum( [ abs(a) for a in diffs ] )
        if sumd > 0:
            print nes > 0
            print "orig", l1
            print "modi", l2
        """

        try:
            top = allHigherAndPos/float(allPos) #p value
            down = nesHigherAndPos/float(nesPos)

            fdrs.append(top/down)
        except:
            fdrs.append(1000000000.0)
    
    return zip(enrichmentScores, nEnrichmentScores, enrichmentPVals, fdrs)

import obiGeneMatch

def nth(l,n): return [ a[n] for a in l ]

def itOrFirst(data):
    if iset(data):
        return data
    else:
        return data[0]

class GSEA(object):

    def __init__(self, organism="hsa"):
        self.a = "blabla"
        self.genesets = {}
        self.organism = organism

    def keepOnlyMeanAttrs(self, datai, atLeast=3, classValues=None):
        """
        Attributes need to be continuous, they need to have
        at least one value.

        In order of attribute to be valid, it needs to have at least
        [atLeast] values for every class value.

        Keep only specified classes - group them in two values.
        """

        cv = itOrFirst(datai).domain.classVar
        nclassvalues = None

        if cv:
            oldcvals = [ a for a in cv.values ]
            
            if not classValues:
                classValues = [ oldcvals[0], oldcvals[1] ]

            toJoin = []

            for vals in classValues:
                if issequencens(vals):
                    toJoin.append(list(vals))
                else:
                    toJoin.append([vals])

            classValues = reduce(lambda x,y: x+y, toJoin)
            classValues = [ str(a) for a in classValues ] # ok class values

            #dictionary of old class -> new class
            mapval = {}
            nclassvalues = [] # need to preserver order

            for joinvals in toJoin:
                joinvalsn = "+".join([ str(val) for val in sorted(joinvals) ])
                nclassvalues.append(joinvalsn)

                for val in joinvals:
                    mapval[str(val)] = joinvalsn

            #take only examples with classValues classes
            nclass = orange.EnumVariable(cv.name, values=nclassvalues)
            ndom = orange.Domain(itOrFirst(datai).domain.attributes, nclass)

            def removeAndTransformClasses(data):
                """
                Removes unnecessary class values and joines them according
                to function input.
                """
                examples = []
                for ex in data:
                    if ex[cv] in classValues:
                        vals = [ ex[a] for a in data.domain.attributes ]
                        vals.append(mapval[str(ex[cv].value)])
                        examples.append(vals)

                return orange.ExampleTable(ndom, examples)

            if iset(datai):
                datai = removeAndTransformClasses(datai)
            else:
                datai = [ removeAndTransformClasses(data) for data in datai ]

            #join specified classes

        def attrOk(a, data):

            #can't
            if a.varType != orange.VarTypes.Continuous:
                return False

            dom2 = orange.Domain([a],data.domain.classVar)
            d2 = orange.ExampleTable(dom2, data)
            vals = [ex[0].value for ex in d2 if not ex[0].isSpecial()]

            if len(vals) < 1:
                return False 
            
            if len(data) > 1 and data.domain.classVar:
                valc = [ [ex[0].value for ex in d2 \
                            if not ex[0].isSpecial() and ex[1] == data.domain.classVar[i] \
                       ] for i in range(len(nclassvalues)) ]
                minl = min( [ len(a) for a in valc ])
                
                if minl < atLeast:
                    return False

            return True
        

        def notOkAttributes(data):
            ignored = []
            for a in data.domain.attributes:
                if not attrOk(a, data):
                    ignored.append(a)
            return ignored
        
        ignored = []
        if iset(datai):
            ignored = set(notOkAttributes(datai))
        else:
            ignored = set(reduce(lambda x,y: x+y, [ notOkAttributes(data) for data in datai ]))

        natts = [ a for a in itOrFirst(datai).domain.attributes if a not in ignored ]
        #print ignored, natts, set(ignored) & set(natts)

        ndom = orange.Domain(natts, itOrFirst(datai).domain.classVar)

        datao = None
        if iset(datai):
            datao = orange.ExampleTable(ndom, datai)
        else:
            datao = [ orange.ExampleTable(ndom, data) for data in datai ]

        return datao, ignored, nclassvalues
        
    def setData(self, data, classValues=None, atLeast=3, caseSensitive=False):

        data, info, classValues  = self.keepOnlyMeanAttrs(data, classValues=classValues, atLeast=atLeast)
        #print "removed attributes", info
        #print "class values taken", classValues

        self.data = data
        attrnames = [ a.name for a in itOrFirst(self.data).domain.attributes ]
        self.gm = obiGeneMatch.GeneMatch(attrnames, organism=self.organism, caseSensitive=caseSensitive)
 
    def addGeneset(self, genesetname, genes):

        #t = time.time()

        if genesetname in self.genesets:
            raise Exception("Geneset with the name " + \
                + genesetname + " is already in genesets.")
        else:
            #print genesetname, genes
            #matching to unified gene names
            datamatch = self.gm.match(genes)
    
            #print "Matched", len(datamatch), "of", len(genes)
            self.genesets[genesetname] = ( genes, datamatch )

        #print "allf", time.time() - t

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
            [(at.name, i) for i,at in enumerate(itOrFirst(data).domain.attributes)])

        #print namesToIndices
        #print self.toRealNames

        nsubsets = []
        nsubsetsNames = []

        for subset in subsetsok:
            nsubsets.append( \
                [namesToIndices[b] for a,b in subset[1][1]])
            nsubsetsNames.append([b for a,b in subset[1][1]])

        if len(nsubsets) == 0:
            return {} # prevent pointless computation of attribe ranks

        #print nsubsets

        if len(itOrFirst(data)) > 1:
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

def getDefaultGenesets():

    def unpckGS(filename):
        import pickle
        f = open(filename,'rb')
        return pickle.load(f)

    import orngRegistry
    return unpckGS(orngRegistry.directoryNames["bufferDir"] + "/gsea/geneSets_MSIGDB.pck")

def runGSEA(data, classValues=None, organism="hsa", geneSets=None, n=100, permutation="class", minSize=3, maxSize=1000, minPart=0.1, atLeast=3, **kwargs):

    gso = GSEA(organism=organism)
    gso.setData(data, classValues=classValues, atLeast=atLeast)
    
    if geneSets == None:
        geneSets = getDefaultGenesets()

    for name,genes in geneSets.items():
        gso.addGeneset(name, genes)

    res1 = gso.compute(n=n, permutation=permutation, minSize=minSize, maxSize=maxSize, minPart=minPart, **kwargs)
    return res1

"""
Genesets
"""

def strip(s):
    #return s
    return s.rstrip().lstrip()

def possiblyReadFile(s):
    """
    If s is not a string, then it is probably a file.
    Read it's contents.
    """
    if isinstance(s, basestring):
        return s
    else:
        return s.read()

def handleNELines(s, fn):
    """
    Run function on nonempty lines of a string.
    Return a list of results for each line.
    """
    lines = s.split("\n")
    lines = [ strip(l) for l in lines ]
    lines = filter(lambda x: x != "", lines)
    return [ fn(l) for l in lines ]

def genesetsLoadGMT(s):
    """
    Eech line consists of tab separated elements. First is
    the geneset name, next is it's description. 
    
    For now the description is skipped.
    """

    s = possiblyReadFile(s)

    def hline(s):
        tabs = [ strip(tab) for tab in s.split("\t") ]
        return tabs[0], tabs[2:]

    return dict(handleNELines(s, hline))

def collectionsPathname():
    if not collectionsPath:
        import orngRegistry
        return os.path.join(orngRegistry.directoryNames["bufferDir"], "gsea", "genesets")
    else:
        return collectionsPath


def createCollection(lnf):
    """
    Input - list of tuples of geneset collection name and GMT
    file.
    """

    def addSource(dic, addition):
        return dict( \
            [ (addition + name,genes) for name,genes in dic.items() ] )

    gen1 = {}

    for n,fn in lnf:
        if fn.lower()[-4:] == ".gmt":
            gen2 = genesetsLoadGMT(open(fn,"rt"))
        elif fn.lower()[-4:] == ".pck":
            import pickle
            f = open(fn,'rb')
            gen2 = pickle.load(f)

        gen1.update(addSource(gen2, "[%s] " % n))

    return gen1

def getCollectionFiles(path=collectionsPathname()):

    def loadInfo(path):
        #TODO load info file
        return {}
        
    info = loadInfo(path)

    def okFile(fn):
        if fn.lower()[-4:] == ".gmt":
            return True
        elif fn.lower()[-4:] == ".pck":
            return True
        return False

    files = sorted(filter(okFile, os.listdir(path)))

    #print files, os.listdir(path)

    out = []

    for file in files:
        fn = os.path.join(path, file)
        name = info.get(file, file)
        if name == file:
            #remove suffix if no info
            if name.lower()[-4:] == ".gmt" or name.lower()[-4:] == ".pck":
                name = name[:-4]

        out.append( (name, fn) )

    return out

def collections(l=[], default=True, path=collectionsPathname()):
    """
    Input is a list of collections.
    Default - if default collections are included.
    Input is a list of names. If names match to any names in path,
    they are taken. If not, file with that name is regarded as
    a filename of gene set colections
    """
    collections = getCollectionFiles(path)

    coln = nth(collections, 0)
    colff = nth(collections, 1)
    colf =  [ os.path.split(f)[1] for f in colff ]

    check = [ coln, colff, colf ]

    choosen = set()
    if default:
        choosen = choosen | set(collections)

    for col in l:
        added = False
        if not added:
            try: # if integer it can be the index
                choosen = choosen | set( [ collections[int(col)] ])
                added = True
            except:
                pass
        if not added:
            for cl in check:
                if col in cl:
                    choosen = choosen | set( [ collections[cl.index(col)] ])
                    added = True
                    break
        if not added:
            choosen = choosen | set( [ (col, col) ] )            

    return createCollection(list(choosen))

"""
End genesets
"""



def etForAttribute(datal,a):
    tables = len(datal)

    def getAttrVals(data, attr):
        dom2 = orange.Domain([data.domain[attr]], False)
        dataa = orange.ExampleTable(dom2, data)
        return [ a[0].native() for a in dataa ]

    domainl = []
    valuesl = []

    for id, data in enumerate(datal):
        v = getAttrVals(data,a)
        valuesl.append(v)
        domainl.append(orange.FloatVariable(name=("v"+str(id))))

    classvals = getAttrVals(data, datal[0].domain.classVar)
    valuesl += [ classvals ]

    dom = orange.Domain(domainl, datal[0].domain.classVar)
    examples = [ list(a) for a in zip(*valuesl) ]

    datat = orange.ExampleTable(dom, examples)

    return datat


def evaluateEtWith(fn):
    """
    fn - evaluates example talbe
    """

    def newf(datal):
        res = []
        for a in datal[0].domain.attributes:
            et = etForAttribute(datal, a)
            res.append(fn(et))
        return res

    return newf


if  __name__=="__main__":

    #import mOrngData

    #data = orange.ExampleTable("iris.tab")
    #data = orange.ExampleTable("DLBCL.tab")
    #data = mOrngData.getAttributes(data, range(100))

    #subsets = [[2,3], [1,4,5,6], [3,4,5], [7,8,9], range(10,20), [14,14,20], [45,51,64], [33,46,49], [31,23,66], [3,54,55], [43,12,96], [1,34,12], [1,43,21], [5,34,87]]

    #gseal = gseaE(data, subsets, rankingFromOrangeMeas(MA_signalToNoise()), n=10, permutation="class")

    #b =  rankingFromOrangeMeas(MA_signalToNoise())(data)
    #gseal = gseaR(b, subsets, n=1000)

    #for el in gseal:
    #    print el
 
    #print gseal

    #collectionsPathname()

    #print collections(["c2.cp.v2.5.symbols.gmt"], default=False)

    #import sys; sys.exit(0)
     
    def unpckGS(filename):
        import pickle
        f = open(filename,'rb')
        return pickle.load(f)

    genesetFile = "genesets/geneSets_cp.pck"
    gen1 = unpckGS(genesetFile)
    
    #data = orange.ExampleTable("leukemiaGSEA.tab")
    #data = orange.ExampleTable("demo.tab")
    data = orange.ExampleTable("sterolTalkHepa.tab")

    print data.domain.classVar.values

    print "loaded data"

    def novi():
        pass
        #print "done"

    def give1(data):
        return 1

    import orngTest, orngStat

    from math import sqrt

    def knn(data):
        learner = orange.kNNLearner(k=int(sqrt(len(data))))
        res = orngTest.crossValidation([ learner ], data)
        return orngStat.AUC(res)[0]**3


    gen1 = collections(["c2.all.v2.5.symbols.gmt"], default=False)
    print gen1.keys()[:100]

    import sys; sys.exit(0)

    #data = orange.ExampleTable(data.domain, data[:1])

    #res2 = runGSEA(data, n=100, geneSets=gen1, permutation="class", callback=novi, atLeast=3)
    res2 = runGSEA([data, data], n=100, geneSets=gen1, permutation="class", callback=novi, atLeast=3, rankingf=evaluateEtWith(knn))
    
    print '\n'.join([ str(a) + ": " +str(b[:3]) for a,b in sorted(res2.items(), key=lambda x: x[1][2])])




