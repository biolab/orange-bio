from __future__ import absolute_import

from collections import defaultdict
import random
import time

import numpy

import orange
import Orange

from . import geneset as obiGeneSets
from .utils.expression import *
from . import gene as obiGene

"""
Gene set enrichment analysis.

Author: Marko Toplak
"""

def iset(data):
    """
    Is data orange.ExampleTable?
    """
    return isinstance(data, orange.ExampleTable)

def issequencens(x):
    "Is x a sequence and not string ? We say it is if it has a __getitem__ method and is not string."
    return hasattr(x, '__getitem__') and not isinstance(x, basestring)

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
    ordered.sort(key=lambda x: -x[1]) #sort by correlation, descending
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

    cors = [ abs(lcor[i])**p for i in subset ] #belowe in numpy
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
        
    last = 0

    maxSum = minSum = csum = 0.0

    for a,b in sorted(map.items()):
        diff = a-last
        csum += notInA*diff
        last = a+1
        
        if csum < minSum:
            minSum = csum
        
        csum += b

        if csum > maxSum:
            maxSum = csum

    #finish it
    diff = (len(ordered))-last
    csum += notInA*diff

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
        n=100, permutation="class", callback=None):
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

    runOptCallbacks(callback)

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

        runOptCallbacks(callback)

    return gseaSignificance(enrichmentScores, enrichmentNulls)


def runOptCallbacks(callback):
    if callback is not None:
        try:
            [ a() for a in callback ]
        except:
            callback()            

def gseaR(rankings, subsets, n, callback=None):
    """
    """
    enrichmentScores = []
    ordered = orderedPointersCorr(rankings)
    
    def rev(l):
        return numpy.argsort(l)

    rev2 = rev(ordered)

    for subset in subsets:

        es = enrichmentScoreRanked(subset, rankings, ordered, rev2=rev2)[0]
        enrichmentScores.append(es)
    
    runOptCallbacks(callback)

    enrichmentNulls = [ [] for a in range(len(subsets)) ]

    for i in range(n):
        
        r2 = shuffleList(rankings, random.Random(2000+i))
        ordered2 = orderedPointersCorr(r2)
        rev22 = rev(ordered2)

        for si,subset in enumerate(subsets):

            esn = enrichmentScoreRanked(subset, r2, ordered2, rev2=rev22)[0]
            enrichmentNulls[si].append(esn)

        runOptCallbacks(callback)

    return gseaSignificance(enrichmentScores, enrichmentNulls)


def gseaSignificance(enrichmentScores, enrichmentNulls):

    #print enrichmentScores

    import time

    tb1 = time.time()

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
 

    #print "First part", time.time() - tb1

    #FDR computation
    #create a histogram of all NES(S,pi) over all S and pi
    vals = reduce(lambda x,y: x+y, nEnrichmentNulls, [])


    def shorten(l, p=10000):
        """
        Take each len(l)/p element, if len(l)/p >= 2.
        """
        e = len(l)/p
        if e <= 1:
            return l
        else:
            return [ l[i] for i in xrange(0, len(l), e) ]

    #vals = shorten(vals) -> this can speed up second part. is it relevant TODO?

    """
    Use this null distribution to compute an FDR q value, for a given NES(S) =
    NES* >= 0. The FDR is the ratio of the percantage of all (S,pi) with
    NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
    observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S)
    = NES* <= 0.
    """

    nvals = numpy.array(sorted(vals))
    nnes = numpy.array(sorted(nEnrichmentScores))

    #print "LEN VALS", len(vals), len(nEnrichmentScores)

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

        #this could be speed up twice with the same accuracy! 
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
    
    #print "Whole part", time.time() - tb1

    return zip(enrichmentScores, nEnrichmentScores, enrichmentPVals, fdrs)


def nth(l,n): return [ a[n] for a in l ]

def itOrFirst(data):
    """ Returns input if input is of type ExampleTable, else returns first
    element of the input list """
    if iset(data):
        return data
    else:
        return data[0]

def wrap_in_list(data):
    """ Wraps orange.ExampleTable in a list """
    if iset(data):
        return [ data ]
    else:
        return data

class transform_class(object):
    
    def __init__(self, cv, mapval, class_values, nclass):
        self.cv = cv
        self.mapval = mapval
        self.class_values = class_values
        self.nclass = nclass
    
    def __call__(self, ex, *args, **kwargs):
        """
        Removes unnecessary class values and joins them according
        to function input.
        """
        if ex[self.cv] in self.class_values:
            return self.nclass(self.mapval[str(ex[self.cv].value)])
        else:
            return "?"

def takeClasses(datai, classValues=None):
    """
    Function joins class groups specified in an input pair
    classValues. Each element of the pair is a list of class
    values to be joined to first or second class. Group
    classes in two new class values.
    If class values are not specified, take all the classes.

    Input data can be a single data set or a list of data sets
    with the same domain.

    Returns transformed data sets / data sets. 
    """

    cv = itOrFirst(datai).domain.classVar
    nclassvalues = None

    if cv and len(itOrFirst(datai)) > 1:
        oldcvals = [ a for a in cv.values ]
        
        if not classValues:
            classValues = oldcvals

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

        tco = transform_class(cv=cv,mapval=mapval,class_values=classValues, nclass=nclass)

        nclass.get_value_from = tco

        ndom = orange.Domain(itOrFirst(datai).domain.attributes, nclass)

        def removeAndTransformClasses(data):
            ndata = Orange.data.Table(ndom, data)
            ndata = Orange.data.Table(ndom, [ex for ex in ndata if ex[-1].value != "?"])
            return ndata

        if iset(datai):
            datai = removeAndTransformClasses(datai)
        else:
            datai = [ removeAndTransformClasses(data) for data in datai ]

    return datai

def removeBadAttributes(datai, atLeast=3):
    """
    Removes attributes which would obscure GSEA analysis.

    Attributes need to be continuous, they need to have
    at least one value. Remove other attributes.

    For the attribute to be valid, it needs to have at least
    [atLeast] values for every class value.

    Return transformed data set / data sets and ignored attributes.
    """

    def attrOk(a, data):
        """
        Attribute is ok if it is continouous and if containg
        at least atLest not unknown values.
        """

        a = data.domain.attributes.index(a)

        #can't
        if data.domain.attributes[a].varType != orange.VarTypes.Continuous:
            return False

        if len(data) == 1:

            vals = [ex[a].value for ex in data if not ex[a].isSpecial()]
            if len(vals) < 1:
                return False 
        
        if len(data) > 1 and data.domain.classVar and atLeast > 0:
            
            dc = dict( (v, 0) for v in data.domain.classVar.values )

            for ex in data:
                if not ex[a].isSpecial():
                    dc[ex[-1].value] += 1

            minl = min(dc.values())

            if minl < atLeast:
                #print "Less than atLeast"
                return False

        return True
    

    def notOkAttributes(data):
        ignored = []
        for a in data.domain.attributes:
            if not attrOk(a, data):
                #print "Removing", a
                ignored.append(a)
        return ignored
    
    ignored = []
    if iset(datai):
        ignored = set(notOkAttributes(datai))
    else:
        #ignore any attribute which is has less than atLeast values for each class
        #ignored = set(reduce(lambda x,y: x+y, [ notOkAttributes(data) for data in datai ]))

        #remove any attribute, which is ok in less than half of the dataset
        ignored = []
        for a in itOrFirst(datai).domain.attributes:
            attrOks = sum([ attrOk(a, data) for data in datai ])
            if attrOks < len(datai)/2:
                ignored.append(a)


    natts = [ a for a in itOrFirst(datai).domain.attributes if a not in ignored ]
    #print ignored, natts, set(ignored) & set(natts)

    ndom = orange.Domain(natts, itOrFirst(datai).domain.classVar)

    datao = None
    if iset(datai):
        datao = orange.ExampleTable(ndom, datai)
    else:
        datao = [ orange.ExampleTable(ndom, data) for data in datai ]

    return datao, ignored

def keepOnlyMeanAttrs(datai, atLeast=3, classValues=None):
    """
    Attributes need to be continuous, they need to have
    at least one value.

    In order of attribute to be valid, it needs to have at least
    [atLeast] values for every class value.

    Keep only specified classes - group them in two values.
    """    
    datai = takeClasses(datai, classValues=classValues)
    return removeBadAttributes(datai, atLeast=atLeast)

def data_single_meas_column(data):
    """ 
    Returns true if data seems to be in one column
    (float variables) only. This column should contain 
    the rankings
    """
    columns = [a for a in data.domain] +  [ data.domain.getmeta(a) for a in list(data.domain.getmetas()) ]
    floatvars = [ a for a in columns if a.varType == orange.VarTypes.Continuous ]
    if len(floatvars) == 1:
        return True
    else:
        return False

def transform_data(data, phenVar, geneVar):
    """
    if we have log2ratio in a single value column, transpose the matrix
    i.e. we have a single column with a continous variable. first
    string variable then becomes the gene name

    The goal is to have different phenotypes annotated with a class,
    and names of genes as attribute names.

    If phenVar is False, then we can work, then the input already
    consists of scores of differential expressions

    If we have a single column, transpose it. 
    If phenVar is one of the groups, transpose the matrix.
    """

    def prepare_data(data, phenVar=None, geneVar=None):

        def rorq(a, name):
            """ Group annatation or question mark. """
            try: 
                return a.attributes[name]
            except: 
                return '?'
   
        #use class as phenotype by default, if it is present,
        #if not, do not use any phenotype!
        if phenVar == None: 
            if not data.domain.classVar:
                phenVar = False
            else:
                phenVar = data.domain.classVar


        #TODO validate phenVar and geneVar?
        #TODO autodetection of groups?

        #transpose is not needed if phenVar is classVar or phenVar is False
        #and there is only one sample
        if isinstance(phenVar, Orange.feature.Descriptor) or \
            (phenVar == False and len(data) == 1):

            if geneVar == None: #if not specified, set as true in this stage
                geneVar = True

            floatvars = [ a for a in data.domain.attributes \
                if a.varType == orange.VarTypes.Continuous ]

            #rename attributes without touching the original variable
            if geneVar != True:
                fl2 = []

                for a in floatvars:
                    na = orange.FloatVariable(name=rorq(a, geneVar))
                    na.getValueFrom = lambda e, rw: e[a]
                    fl2.append(na)

                floatvars = fl2

            dom = orange.Domain(floatvars, phenVar)
            return orange.ExampleTable(dom, data)

        elif phenVar == False or isinstance(phenVar, basestring):

            cands = allgroups(data)
            pv = False
            if phenVar != False:
                pv = orange.EnumVariable(name="phenotype", 
                    values=list(cands[phenVar]))

            #take the only string attribute as a gene name
            gc = gene_cands(data, False)
            if geneVar == None:
                if len(gc) == 1:
                    geneVar = gc[0]
                else:
                    geneNamesUnspecifiedError()
           
            latts = [ orange.FloatVariable(name=ex[geneVar].value) \
                for ex in data ]

            domain = orange.Domain(latts, pv)

            examples = []
            for at in data.domain.attributes:
                if at.varType == orange.VarTypes.Continuous:
                    vals = [ ex[at].value for ex in data ]
                    if pv != False: #add class value
                        vals.append(rorq(at, phenVar))
                    examples.append(orange.Example(domain, vals))

            return orange.ExampleTable(domain, examples)
        else:
            wrongInputsError()

    #transform all example tables
    single = iset(data)
    transposed = [ prepare_data(d, phenVar, geneVar) for d in wrap_in_list(data) ]

    if single:
        return transposed[0]
    else:
        return transposed


def allgroups(data):
    """
    Return all phenotype descriptors of attributes with their values.
    """
    sd = defaultdict(set)
    for attr in data.domain.attributes:
        for key, value in attr.attributes.items():
            sd[key].add(value)
    return sd

def already_have_correlations(data):
    return len(data) == 1 or data_single_meas_column(data)

def need_to_transpose_single(data):
    if len(data) == 1:
        return False
    else:
        return True

def phenotype_cands(data):
    """
    Return all phenotype candidate descriptors in a list of tuples 
    (variable, values). Candidates are class variable, if it exists and 
    attributes dictionaries of attributes.
    Phenotype candidates must contain at least two differend values.
    """
    if already_have_correlations(data):
        return [ (False, set()) ]
    else:
        cv = []
        if data.domain.classVar and data.domain.classVar.varType == orange.VarTypes.Discrete:
            cv.append((data.domain.classVar, set(data.domain.classVar.values)))
        cands = cv + sorted(allgroups(data).items())
        return filter(lambda x: len(x[1]) >= 2, cands)

def gene_cands(data, correct):
    """
    Returns all valid gene descriptors with regards to the choosen
    phenotype variable.
    Return variable descriptor for variables, name of the group for
    descriptions in attr.attributes and True for the usage
    of attribute names.
    Correct is True, if the example table has genes as attributes.
    """
    if correct:
        #gene names could be in attributes or as gene names (marker True)
        return [True] + nth(sorted(allgroups(data)),0)
    else:
        #gene names are values of some string attribute
        columns = [a for a in data.domain] +  \
            [ data.domain.getmeta(a) for a in list(data.domain.getmetas()) ]
        stringvars = [ a for a in columns if a.varType == 6 ]
        return stringvars

def is_variable(phenVar):
    return isinstance(phenVar, orange.Variable)

class GSEA(object):

    def __init__(self, data, organism=None, matcher=None, classValues=None, 
        atLeast=3, caseSensitive=False, phenVar=None, geneVar=None):
        """
        If the data set constains multiple measurements for a single gene,
        all are considered. Individual constributions of such measurements
        are not weighted down - each measurement is as important as they
        would measure different genes.

        phenVar and geneVar can ether be an orange attribute or a string.
        If they are strings, then they describe a group.
        """

        self.genesets = {}
        self.organism = organism

        if organism != None:
            print "WARNING: obiGsea - organism and caseSensitive parameters are deprecated. Use matcher instead."

        self.gsweights = {}
        self.namesToIndices = None
        self.gm = matcher

        data = transform_data(data, phenVar, geneVar)

        data, info = keepOnlyMeanAttrs(data, classValues=classValues, atLeast=atLeast)

        self.data = data

        #init attrnames
        attrnames = [ a.name for a in itOrFirst(self.data).domain.attributes ]

        if self.gm == None: #build a gene matcher, if if does not exists
            self.gm = obiGene.matcher([obiGene.GMKEGG(self.organism, ignore_case=not caseSensitive)], 
                ignore_case=not caseSensitive, direct=True)
            print "WARNING: gene matcher build automatically for organism: " + self.organism

        self.gm.set_targets(attrnames)

 
    def addGeneset(self, genesetname, genes):
        """
        Add a single gene set. See addGenesets function.
        Solely for backwards compatibility.
        """
        self.addGenesets({ genesetname: genes })

    def addGenesets(self, genesets):
        """
        Adds genesets from input dictionary. Performs gene matching. Adds
        to a self.genesets: key is genesetname, it's values are individual
        genes and match results.
        """
        for g in obiGeneSets.GeneSets(genesets):
            genes = g.genes
            datamatch = filter(lambda x: x[1] != None, 
                [ (gene, self.gm.umatch(gene)) for gene in genes])
            self.genesets[g] = datamatch

    def selectGenesets(self, minSize=3, maxSize=1000, minPart=0.1):
        """ Returns a list of gene sets that have sizes in limits """

        def okSizes(orig, transl):
            """compares sizes of genesets to limitations"""
            if len(transl) >= minSize and len(transl) <= maxSize \
                and float(len(transl))/len(orig) >= minPart:
                return True
            return False

        return  dict( (a,c) for a,c in self.genesets.iteritems() if okSizes(a.genes,c) )

    def genesIndices(self, genes):
        """
        Returns in attribute indices of given genes.
        Buffers locations dictionary.
        """
        if not self.namesToIndices:
            self.namesToIndices = defaultdict(list)
            for i,at in enumerate(itOrFirst(self.data).domain.attributes):
                self.namesToIndices[at.name].append(i)
        return reduce(lambda x,y:x+y, [ self.namesToIndices[gname] for gname in genes ], [])

    def to_gsetsnum(self, gsets):
        """
        Returns a dictionary of given  gene sets in gsetnums format.
        """
        return dict( (gs, self.genesIndices(nth(self.genesets[gs],1))) for gs in gsets)

    def compute(self, minSize=3, maxSize=1000, minPart=0.1, n=100, callback=None, rankingf=None, permutation="class"):

        subsetsok = self.selectGenesets(minSize=minSize, maxSize=maxSize, minPart=minPart)

        gsetsnum = self.to_gsetsnum(subsetsok.keys())
        gsetsnumit = gsetsnum.items() #to fix order

        #gsetsnumit = gsetsnumit[:1]
        #print gsetsnumit

        if len(gsetsnum) == 0:
            return {} # quick return if no genesets

        if len(itOrFirst(self.data)) > 1:
            gseal = gseaE(self.data, nth(gsetsnumit,1), n=n, callback=callback, permutation=permutation, rankingf=rankingf)
        else:
            rankings = [ self.data[0][at].native() for at in self.data.domain.attributes ]
            gseal = gseaR(rankings, nth(gsetsnumit,1), n, callback=None)

        res = {}

        for gs, gseale in zip(nth(gsetsnumit,0),gseal):
            rdict = {}
            rdict['es'] = gseale[0]
            rdict['nes'] = gseale[1]
            rdict['p'] = gseale[2]
            rdict['fdr'] = gseale[3]
            rdict['size'] = len(gs.genes)
            rdict['matched_size'] = len(self.genesets[gs])
            rdict['genes'] = nth(self.genesets[gs],1)
            res[gs] = rdict

        return res

def direct(data, gene_sets, matcher, min_size=3, max_size=1000, min_part=0.1,
    gene_desc=None, n=100, callback=None):
    """ Gene Set Enrichment analysis for pre-computed correlations
    between genes and phenotypes. 
    
    :param Orange.data.Table: Precomputed correlations as a data set 
        with a single continuous variable or a single 
        :obj:`~Orange.data.Instance`.

    See :obj:`run` for other parameters.

    :return: See :obj:`run`. 
    """

    assert len(data.domain.attributes) == 1 or len(data) == 1
    return runGSEA(data, geneSets=gene_sets, matcher=matcher, minSize=min_size, 
        maxSize=max_size, minPart=min_part, n=n, geneVar=gene_desc, callback=callback)

def run(data, gene_sets, matcher, min_size=3, max_size=1000, min_part=0.1,
    at_least=3, phenotypes=None, gene_desc=None, phen_desc=None, n=100, 
    permutation="phenotype", callback=None, rankingf=None):
    """ Run Gene Set Enrichment Analysis.

    :param Orange.data.Table data: Gene expression data.  
    :param Orange.bio.geneset.GeneSets gene_sets: Gene sets.  
    :param Orange.bio.gene.Matcher matcher: Initialized gene matcher.
    :param tuple phenotypes: A pair describing two distinct phenotypes.
        Each element can also be a list of values. Only examples with
        a chosen phenotypes are analysed. Default: values phenotypes
        of ``phen_desc``.
    :param n: Number of permutations for significance computation. Default: 100.
    :param str permutation: Permutation type, "phenotype" (default) for 
        phenotypes, "gene" for genes.
    :param int min_size:
    :param int max_size: Minimum and maximum allowed number of genes from
        gene set also the data set. Defaults: 3 and 1000.
    :param float min_part: Minimum fraction of genes from the gene set
        also in the data set. Default: 0.1.
    :param int at_least: Minimum number of valid gene values for each 
        phenotype (the rest are ignored). Default: 3.
    :param phen_desc: Location of data on phenotypes. By default the
        ``data.domain.class_var`` is used if it exists. If string, the
        corresponding entry from ``attributes`` dictionary of individual
        features specifies the phenotype. In latter case, each attribute
        represents one sample.
    :param gene_desc: Locations of gene names. If True, gene names
        are attribute names. If a string, that entry of the individual
        features'``attributes`` dictionary is used. If each attribute
        specifies a sample, then the user should pass the meta variable
        containing the gene names. Defaults to attribute names if each
        example specifies one sample.

    :return: | a dictionary where key is a gene set and values are:
        | { es: enrichment score, 
        | nes: normalized enrichment score, 
        | p: P-value, 
        | fdr: FDR, 
        | size: gene set size,
        | matched_size: genes matched to the data, 
        | genes: gene names from the data set }

    """
    assert len(data.domain.attributes) > 1 or len(data) > 1
    assert permutation in ["phenotype", "gene"]
    if permutation == "phenotype":
        permutation = "class"
    return runGSEA(data, geneSets=gene_sets, matcher=matcher, minSize=min_size, 
        maxSize=max_size, minPart=min_part, n=n, permutation=permutation, 
        geneVar=gene_desc, callback=callback, phenVar=phen_desc, 
        classValues=phenotypes)

def runGSEA(data, organism=None, classValues=None, geneSets=None, n=100, 
        permutation="class", minSize=3, maxSize=1000, minPart=0.1, atLeast=3, 
        matcher=None, geneVar=None, phenVar=None, caseSensitive=False, 
        rankingf=None, callback=None):
    gso = GSEA(data, organism=organism, matcher=matcher, 
        classValues=classValues, atLeast=atLeast, caseSensitive=caseSensitive,
        geneVar=geneVar, phenVar=phenVar)
    gso.addGenesets(geneSets)
    res1 = gso.compute(n=n, permutation=permutation, minSize=minSize,
        maxSize=maxSize, minPart=minPart, rankingf=rankingf,
        callback=callback)
    return res1

def etForAttribute(datal,a):
    """
    Builds an example table for a single attribute across multiple 
    example tables.
    """

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


def evaluateEtWith(fn, *args, **kwargs):
    """
    fn - evaluates example table given
    following arguments.
    """

    def newf(datal):
        res = []
        for a in datal[0].domain.attributes:
            et = etForAttribute(datal, a)
            res.append(fn(et, *args, **kwargs))
        return res

    return newf


def hierarchyOutput(results, limitGenes=50):
    """
    Transforms results for use by hierarchy output from GO.

    limitGenes - maximum number of genes on output.
    """
    trans = []
    
    for name, res in results.items():
        try:
            second = name.split(' ')[2]
            name = second if second[:2] == 'GO' else name
        except:
            pass
        
        trans.append((name, abs(res["nes"]), res["matched_size"], res["size"], res["p"], min(res["fdr"], 1.0), res["genes"][:limitGenes]))

    return trans

if  __name__=="__main__":

    """
    Old example with another measure function
    data = orange.ExampleTable("gene_three_lines_log.tab")
    data = orange.ExampleTable("sterolTalkHepa.tab")
    gen1 = collections(['steroltalk.gmt', ':kegg:hsa'], default=False)
    rankingf = rankingFromOrangeMeas(MA_anova())
    matcher = obiGene.matcher([obiGene.GMKEGG('hsa')])
    geneVar = gene_cands(data, False)[1]
    phenVar = "group"
    out = runGSEA(data, n=10, geneSets=gen1, permutation="gene", atLeast=3, matcher=matcher, rankingf=rankingf, phenVar=phenVar, geneVar=geneVar)
    """
    """
    data = orange.ExampleTable("sterolTalkHepa.tab")
    gen1 = obiGeneSets.collections('steroltalk.gmt', (("KEGG",), "9606"))
    matcher = obiGene.matcher([obiGene.GMKEGG('hsa')])
    out = runGSEA(data, n=10, geneSets=gen1, permutation="gene", atLeast=3, matcher=matcher)
    """
    matcher = obiGene.matcher([[obiGene.GMKEGG('ddi'), obiGene.GMDicty()]])
    data = orange.ExampleTable("/home/marko/t_gsea1.tab")
    gen1 = obiGeneSets.collections((("KEGG",), "352472"))
    out = runGSEA(data, n=10, geneSets=gen1, permutation="gene", atLeast=3, matcher=matcher, phenVar="growth")
    print out
    print "\n".join(map(str,sorted(out.items())))
    
