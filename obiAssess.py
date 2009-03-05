import obiGsea
import obiGeneSets
import orange
import stats
import statc
import numpy
import math

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
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, rankingf=None):
       
        gso = obiGsea.GSEA(organism=organism)
        gso.setData(data, classValues=classValues, atLeast=0)

        for name,genes in geneSets.items():
            gso.addGeneset(name, genes)

        if rankingf == None:
            rankingf = AT_edelmanParametricLearner()

        oknames = gso.selectGenesets(minSize=minSize, maxSize=maxSize, minPart=minPart).keys()
        gsetsnum = gso.to_gsetsnum(oknames)

        attrans = [ rankingf(iat, gso.data) for iat, at in enumerate(data.domain.attributes) ]

        return Assess(attrans=attrans, gsetsnum=gsetsnum, genesets=gso.genesets)

if __name__ == "__main__":
    
    data = orange.ExampleTable("sterolTalkHepa.tab")
    a = AssessLearner()
    ass = a(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False), rankingf=AT_loessLearner())

    ar = {}

    for d in data:
        print data.domain.classVar.values
        if d[-1].value in ["LK935_48h", "Rif_12h"]:
            cr = ass(d)
            for a,b in cr.items():
                l = ar.get(a, [])
                l.append(b)
                ar[a] = l

    ol =  sorted(ar.items())

    print '\n'.join([ str(a) + ": " +str(b) for a,b in ol])
