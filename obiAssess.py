"""
Construction of gene set scores for each sample.

Learners in this module build models needed for construction
of features from individual genes.

The other classes just take example and return a
dictionary of { name: score } for that example.
"""

import obiGsea
import obiGeneSets
import orange
import stats
import statc
import numpy
import math
import obiExpression
import orngRegression

def normcdf(x, mi, st):
    return 0.5*(2. - stats.erfcc((x - mi)/(st*math.sqrt(2))))

class AT_edelmanParametric(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, nval):

        if self.mi1 == None or self.mi2 == None or self.st1 == None or self.st2 == None:
            return 0 

        try:
            val = nval.value
            if nval.isSpecial():
                return 0 
        except:
            val = nval

        try:
            if val >= self.mi1:
                p1 = 1 - normcdf(val, self.mi1, self.st1)
            else:
                p1 = normcdf(val, self.mi1, self.st1)

            if val >= self.mi2:
                p2 = 1 - normcdf(val, self.mi2, self.st2)
            else:
                p2 = normcdf(val, self.mi2, self.st2)

            #print p1, p2
            return math.log(p1/p2)
        except:
            #print p1, p2, "exception"
            return 0

class AT_edelmanParametricLearner(object):
    """
    Returns attribute transfromer for Edelman parametric measure for a given attribute in the
    dataset.
    Edelman et al, 06.

    Modified a bit.
    """

    def __init__(self, a=None, b=None):
        """
        a and b are choosen class values.
        """
        self.a = a
        self.b = b

    def __call__(self, i, data):
        cv = data.domain.classVar
        #print data.domain

        if self.a == None: self.a = cv.values[0]
        if self.b == None: self.b = cv.values[1]

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[-1].value == value and not ex[i].isSpecial() ]

        list1 = avWCVal(self.a)
        list2 = avWCVal(self.b)

        mi1 = mi2 = st1 = st2 = None

        try:
            mi1 = statc.mean(list1)
            st1 = statc.std(list1)
        except:
            pass
    
        try:
            mi2 = statc.mean(list2)
            st2 = statc.std(list2)
        except:
            pass

        return AT_edelmanParametric(mi1=mi1, mi2=mi2, st1=st1, st2=st2)

class AT_loess(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, nval):

        val = nval.value
        if nval.isSpecial():
            return 0.0 #middle value
        #return first class probablity

        import math

        def saveplog(a,b):
            try:
                return math.log(a/b)
            except:
                if a < b:
                    return -10
                else:
                    return +10

        try:
            ocene = self.condprob(val)
            if sum(ocene) < 0.01:
                return 0.0
            return saveplog(ocene[0], ocene[1])

        except:
            return 0.0

class AT_loessLearner(object):

    def __call__(self, i, data):
        cv = data.domain.classVar
        #print data.domain
        try:
            ca = orange.ContingencyAttrClass(data.domain.attributes[i], data)
            a = orange.ConditionalProbabilityEstimatorConstructor_loess(ca, nPoints=5)
            return AT_loess(condprob=a)
        except:
            return AT_loess(condprob=None)

def nth(l, n):
    return [a[n] for a in l]

class Assess(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        enrichmentScores = [] 

        lcor = [ self.attrans[at](example[at]) for at in range(len(self.attrans)) ]

        ordered = obiGsea.orderedPointersCorr(lcor)

        def rev(l):
           return numpy.argsort(l)

        rev2 = rev(ordered)

        gsetsnumit = self.gsetsnum.items()

        enrichmentScores = dict( 
            (name, obiGsea.enrichmentScoreRanked(subset, lcor, ordered, rev2=rev2)[0]) \
            for name,subset in gsetsnumit)
    
        return enrichmentScores


class AssessLearner(object):
    """
    Uses the underlying GSEA code to select genes.
    Takes data and creates attribute transformations.
    """
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, rankingf=None):
        data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        
        if rankingf == None:
            rankingf = AT_edelmanParametricLearner()

        #rank individual attributes on the training set
        attrans = [ rankingf(iat, data) for iat, at in enumerate(data.domain.attributes) ]

        return Assess(attrans=attrans, gsetsnum=gsetsnum)

def selectGenesetsData(data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
    """
    Returns gene sets and data which falling under upper criteria.
    """
    gso = obiGsea.GSEA(data, organism=organism, classValues=classValues, atLeast=0)
    gso.addGenesets(geneSets)
    oknames = gso.selectGenesets(minSize=minSize, maxSize=maxSize, minPart=minPart).keys()
    gsetsnum = gso.to_gsetsnum(oknames)
    return gso.data, oknames, gsetsnum

def ideker_activity_score(ex, corg):
    """ activity score for a sample for pathway given by corgs """
    #print [ ex[i].value for i in corg ] #FIXME what to do with unknown values?
    return sum(ex[i].value if ex[i].value != '?' else 0.0 for i in corg)/len(corg)**0.5

class Ideker(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        return dict( (name,ideker_activity_score(example, corg)) \
            for name, corg in self.corgs.items() )

class IdekerLearner(object):
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
        """
        WARNING: input has to be z_ij table! each gene needs to be normalized
        (mean=0, stdev=1) for all samples.
        """
        data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
    
        tscorecache = {}

        def tscorec(data, at, cache=None):
            """ Cached attribute  tscore calculation """
            if cache != None and at in cache: return cache[at]
            ma = obiExpression.MA_t_test()(at,data)
            if cache != None: cache[at] = ma
            return ma

        def compute_corg(data, inds):
            """
            Compute CORG for this geneset specified with gene inds
            in the example table. Output is the list of gene inds
            in CORG.

            """
            #order member genes by their t-scores: decreasing, if av(t-score) >= 0,
            #else increasing
            tscores = [ tscorec(data, at, tscorecache) for at in inds ]
            sortedinds = nth(sorted(zip(inds,tscores), key=lambda x: x[1], \
                reverse=statc.mean(tscores) >= 0), 0)

            def S(corg):
                """ Activity score separation - S(G) in 
                the article """
                asv = orange.FloatVariable(name='AS')
                asv.getValueFrom = lambda ex,rw: orange.Value(asv, ideker_activity_score(ex, corg))
                data2 = orange.ExampleTable(orange.Domain([asv], data.domain.classVar), data)
                return abs(tscorec(data2, 0)) #FIXME absolute - nothing in the article about it
                    
            #greedily find CORGS procing the beset separation
            g = S(sortedinds[:1])
            bg = 1
            for a in range(2, len(sortedinds)+1):
                tg = S(sortedinds[:a])
                if tg > g:
                    g = tg
                    bg = a
                else:
                    break
                
            return sortedinds[:a]

        #now, on the learning set produce the CORGS for each geneset and
        #save it for use in further prediction

        corgs = {}

        for name, inds in gsetsnum.items():
            inds = sorted(set(inds)) # take each gene only once!
            corgs[name] = compute_corg(data, inds)

        return Ideker(corgs=corgs)

def pls_transform(example, constt):
    inds, xmean, xstd, W = constt
    dom = orange.Domain([example.domain.attributes[i1] for i1 in inds ], False)
    newex = orange.ExampleTable(dom, [example])
    ex = newex.toNumpy()[0]
    print "ex", ex
    ex = (ex - xmean) / xstd
    print "ex", ex
    print "W", W
    return numpy.dot(ex,W)

class PLS(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        return dict( (name, pls_transform(example, constt)) \
            for name, constt in self.constructt.items() )

class PLSLearner(object):
    """ Transforms gene sets using Principal Leasts Squares. """
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):

        data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
    
        constructt = {}

        #build weight matrices for every gene set
        for name, inds in gsetsnum.items():
            print name, inds
            dom2 = orange.Domain([ data.domain.attributes[i1] for i1 in inds ], data.domain.classVar)
            data_gs = orange.ExampleTable(dom2, data)
            pls = orngRegression.PLSRegressionLearner(data_gs, nc=1, y=[data_gs.domain.classVar], save_partial=True)
            constructt[name] = inds, pls.XMean, pls.XStd, pls.W
            pt1 = pls.T[0]
            print "T", pls.T
            print "T1_1", pt1
            print data_gs[0]
            pt2 = pls_transform(data[0], constructt[name])
            print "T1_2", pt2

        return PLS(constructt=constructt)

class SimpleFun(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        #for  name,ids in self.gsets.items():
        #    print name, [example[i].value for i in ids], self.fn([example[i].value for i in ids])
        return dict( (name, self.fn([example[i].value for i in ids])) \
            for name,ids in self.gsets.items() )

class SimpleFunLearner(object):
    """
    Just applies a function taking attribute values of an example and
    produces to all gene sets.    
    """
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, fn=None):
        data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        return SimpleFun(gsets=gsetsnum, fn=fn)

class MedianLearner(object):
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, fn=None):
       sfl =  SimpleFunLearner()
       return sfl(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues, fn=statc.median)

class MeanLearner(object):
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, fn=None):
       sfl =  SimpleFunLearner()
       return sfl(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues, fn=statc.mean)

def impute_missing(data):
    #remove attributes with only unknown values
    newatts = []
    for at in data.domain.attributes:
        svalues = [ 1 for a in data if a[at].isSpecial() ]
        real = len(data) - len(svalues)
        if real > 0:
            newatts.append(at)

    dom2 = orange.Domain(newatts, data.domain.classVar)
    data = orange.ExampleTable(dom2, data)

    #impute
    import orngTree 
    imputer = orange.ImputerConstructor_model() 
    imputer.learnerContinuous = imputer.learnerDiscrete = orange.MajorityLearner()
    imputer = imputer(data)

    data = imputer(data)
    return data


if __name__ == "__main__":
    
    data = orange.ExampleTable("sterolTalkHepa.tab")
    data = impute_missing(data)

    #ass = AssessLearner()(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False), rankingf=AT_loessLearner())
    #ass = MeanLearner()(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False))
    ass = PLSLearner()(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False), classValues=["LK935_48h", "Rif_12h"])

    fsdfsd()

    ar = {}


    print data.domain.classVar.values

    for d in data:
        if d[-1].value in ["LK935_48h", "Rif_12h"]:
            cr = ass(d)
            for a,b in cr.items():
                l = ar.get(a, [])
                l.append(b)
                ar[a] = l

    ol =  sorted(ar.items())

    #print '\n'.join([ str(a) + ": " +str(b) for a,b in ol])
