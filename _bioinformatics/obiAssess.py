"""
Construction of gene set scores for each sample.

Learners in this module build models needed for construction
of features from individual genes.

The other classes just take example and return a
dictionary of { name: score } for that example.
"""

from __future__ import absolute_import

from collections import defaultdict
import math

import numpy

import orange, Orange, statc

from . import obiExpression, obiGene, obiGsea, obiGeneSets, stats

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
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, rankingf=None):
        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        
        if rankingf == None:
            rankingf = AT_edelmanParametricLearner()

        #rank individual attributes on the training set
        attrans = [ rankingf(iat, data) for iat, at in enumerate(data.domain.attributes) ]

        return Assess(attrans=attrans, gsetsnum=gsetsnum)

def selectGenesetsData(data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
    """
    Returns gene sets and data which falling under upper criteria.
    """
    gso = obiGsea.GSEA(data, matcher=matcher, classValues=classValues, atLeast=0)
    gso.addGenesets(geneSets)
    okgenesets = gso.selectGenesets(minSize=minSize, maxSize=maxSize, minPart=minPart).keys()
    gsetsnum = gso.to_gsetsnum(okgenesets)
    return gso.data, okgenesets, gsetsnum

def corgs_activity_score(ex, corg):
    """ activity score for a sample for pathway given by corgs """
    #print [ ex[i].value for i in corg ] #FIXME what to do with unknown values?
    return sum(ex[i].value if ex[i].value != '?' else 0.0 for i in corg)/len(corg)**0.5

class CORGs(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        return dict( (name,corgs_activity_score(example, corg)) \
            for name, corg in self.corgs.items() )

class CORGsLearner(object):
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
        """
        WARNING: input has to be z_ij table! each gene needs to be normalized
        (mean=0, stdev=1) for all samples.
        """
        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
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
                asv.getValueFrom = lambda ex,rw: orange.Value(asv, corgs_activity_score(ex, corg))
                data2 = orange.ExampleTable(orange.Domain([asv], data.domain.classVar), data)
                return abs(tscorec(data2, 0)) #FIXME absolute - nothing in the article about it
                    
            #greedily find CORGS procing the best separation
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

        return CORGs(corgs=corgs)

class GSA(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        return dict( (name, statc.mean([example[i].value for i in inds if example[i].value != "?"]) ) \
            for name, inds in self.subsets.items() )

class GSALearner(object):
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
        """
        WARNING: input has to be z_ij table! each gene needs to be normalized
        (mean=0, stdev=1) for all samples.
        """
        import scipy.stats

        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
    
        def tscorec(data, at, cache=None):
            ma = obiExpression.MA_t_test()(at,data)
            return ma

        tscores = [ tscorec(data, at) for at in data.domain.attributes ]

        def to_z_score(t):
            return float(scipy.stats.norm.ppf(scipy.stats.t.cdf(t, len(data)-2)))

        zscores = map(to_z_score, tscores)

        subsets = {}

        for name, inds in gsetsnum.items():
            inds = sorted(set(inds)) # take each gene only once!

            D = statc.mean([max(zscores[i],0) for i in inds]) \
                + statc.mean([min(zscores[i],0) for i in inds])

            if D >= 0:
                subsets[name] = [ i for i in inds if zscores[i] > 0.0 ]
            else:
                subsets[name] = [ i for i in inds if zscores[i] < 0.0 ]

        return GSA(subsets=subsets)

def pls_transform(example, constt):
    """
    Uses calculated PLS weights to transform the given example.
    """

    inds, xmean, W, P = constt
    dom = orange.Domain([example.domain.attributes[i1] for i1 in inds ], False)
    newex = orange.ExampleTable(dom, [example])
    ex = newex.toNumpy()[0]
    ex = ex - xmean # same input transformation

    nc = W.shape[1]

    TR = numpy.empty((1, nc))
    XR = ex

    dot = numpy.dot

    for i in range(nc):
       t = dot(XR, W[:,i].T)
       XR = XR - t*numpy.array([P[:,i]])
       TR[:,i] = t

    return TR

def PLSCall(data, y=None, x=None, nc=None, weight=None, save_partial=False):

    def normalize(vector):
        return vector / numpy.linalg.norm(vector)

    if y == None:
        y = [ data.domain.classVar ]
    if x == None:
        x = [v for v in data.domain.variables if v not in y]

    Ncomp = nc if nc is not None else len(x)
        
    dataX = orange.ExampleTable(orange.Domain(x, False), data)
    dataY = orange.ExampleTable(orange.Domain(y, False), data)

    # transformation to numpy arrays
    X = dataX.toNumpy()[0]
    Y = dataY.toNumpy()[0]

    # data dimensions
    n, mx = numpy.shape(X)
    my = numpy.shape(Y)[1]

    # Z-scores of original matrices
    YMean = numpy.mean(Y, axis = 0)
    XMean = numpy.mean(X, axis = 0)
    
    X = (X-XMean)
    Y = (Y-YMean)

    P = numpy.empty((mx,Ncomp))
    T = numpy.empty((n,Ncomp))
    W = numpy.empty((mx,Ncomp))
    E,F = X,Y

    dot = numpy.dot
    norm = numpy.linalg.norm

    #PLS1 - from Gutkin, shamir, Dror: SlimPLS

    for i in range(Ncomp):
        w = dot(E.T,F)
        w = w/norm(w) #normalize w in Gutkin et al the do w*c, where c is 1/norm(w)
        t = dot(E, w) #t_i -> a row vector
        p = dot(E.T, t)/dot(t.T, t) #p_i t.T is a row vector - this is inner(t.T, t.T)
        q = dot(F.T, t)/dot(t.T, t) #q_i
            
        E = E - dot(t, p.T)
        F = F - dot(t, q.T)

        T[:,i] = t.T
        W[:,i] = w.T
        P[:,i] = p.T

    return XMean, W, P, T

class PLS(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):

        od = {}

        for name, constt in self.constructt.items():
            ts = pls_transform(example, constt)[0]
            if len(ts) == 1:
                od[name] = ts[0]
            else:
                for i,s in enumerate(ts):
                   od[name + "_LC_" + str(i+1)] = s
 
        return od

class PLSLearner(object):
    """ Transforms gene sets using Principal Leasts Squares. """
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, components=1):
        """
        If more that 1 components are used, _LC_componetsNumber is appended to 
        the name of the gene set.
        """

        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
    
        constructt = {}

        #build weight matrices for every gene set
        for name, inds in gsetsnum.items():
            dom2 = orange.Domain([ data.domain.attributes[i1] for i1 in inds ], data.domain.classVar)
            data_gs = orange.ExampleTable(dom2, data)
            xmean, W, P, T = PLSCall(data_gs, nc=components, y=[data_gs.domain.classVar], save_partial=True)
            constructt[name] = inds, xmean, W, P

        return PLS(constructt=constructt)

class PCA(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        od = {}

        for name, constt in self.constructt.items():
            ts = pca_transform(example, constt)[0]
            od[name] = ts

        return od

def pca_transform(example, constt):
    inds, evals, evect, xmean = constt
    dom = orange.Domain([example.domain.attributes[i1] for i1 in inds ], False)
    newex = orange.ExampleTable(dom, [example])
    arr = newex.toNumpy()[0]
    arr = arr - xmean # same input transformation - a row in a matrix

    ev0 = evect[0] #this is a row in a matrix - do a dot product
    a = numpy.dot(arr, ev0)
    return a

def pca(data, snapshot=0):
    "Perform PCA on M, return eigenvectors and eigenvalues, sorted."
    M = data.toNumpy("a")[0]
    XMean = numpy.mean(M, axis = 0)
    M = M - XMean

    T, N = numpy.shape(M)
    # if there are less rows T than columns N, use snapshot method
    if (T < N) or snapshot:
        C = numpy.dot(M, numpy.transpose(M))
        evals, evecsC = numpy.linalg.eigh(C) #columns of evecsC are eigenvectors
        evecs = numpy.dot(M.T, evecsC)/numpy.sqrt(numpy.abs(evals))
    else:
        K = numpy.dot(numpy.transpose(M), M)
        evals, evecs = numpy.linalg.eigh(K)
    
    evecs = numpy.transpose(evecs)

    # sort the eigenvalues and eigenvectors, decending order
    order = (numpy.argsort(numpy.abs(evals))[::-1])
    evecs = numpy.take(evecs, order, 0)
    evals = numpy.take(evals, order)
    return evals, evecs, XMean

class PCALearner(object):
    """ Transforms gene sets using Principal Leasts Squares. """
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):

        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
    
        constructt = {}

        #build weight matrices for every gene set
        for name, inds in gsetsnum.items():
            dom2 = orange.Domain([ data.domain.attributes[i1] for i1 in inds ], data.domain.classVar)

            data_gs = orange.ExampleTable(dom2, data)
            evals, evect, xmean = pca(data_gs)
            constructt[name] = inds, evals, evect, xmean

        return PCA(constructt=constructt)


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
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, fn=None):
        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        return SimpleFun(gsets=gsetsnum, fn=fn)

class MedianLearner(object):
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, fn=None):
       sfl =  SimpleFunLearner()
       return sfl(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues, fn=statc.median)

class MeanLearner(object):
    
    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, fn=None):
       sfl =  SimpleFunLearner()
       return sfl(data, matcher, geneSets, \
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
    from Orange.orng import orngTree 
    imputer = orange.ImputerConstructor_model() 
    imputer.learnerContinuous = imputer.learnerDiscrete = orange.MajorityLearner()
    imputer = imputer(data)

    data = imputer(data)
    return data

def setSig_example(ldata, ex, genesets):
    """
    Create an dictionary with (geneset_name, geneset_expression) for every
    geneset on input.

    ldata is class-annotated data
    """
    #seznam ocen genesetov za posamezen primer v ucni mzozici
    pathways = {}

    def setSig_example_geneset(ex, data):
        """ ex contains only selected genes """

        distances = [ [], [] ]    

        def pearsonr(v1, v2):
            try:
                return statc.pearsonr(v1, v2)[0]
            except:
                return numpy.corrcoef([v1, v2])[0,1]

        def pearson(ex1, ex2):
            attrs = range(len(ex1.domain.attributes))
            vals1 = [ ex1[i].value for i in attrs ]
            vals2 = [ ex2[i].value for i in attrs ]
            return pearsonr(vals1, vals2)

        def ttest(ex1, ex2):
            try:
                return stats.lttest_ind(ex1, ex2)[0]
            except:
                return 0.0
        
        #maps class value to its index
        classValueMap = dict( [ (val,i) for i,val in enumerate(data.domain.classVar.values) ])
     
        #create distances to all learning data - save or other class
        for c in data:
            distances[classValueMap[c[-1].value]].append(pearson(c, ex))

        return ttest(distances[0], distances[1])
           
    for name, gs in genesets.items(): #for each geneset
        #for each gene set: take the attribute subset and work on the attribute subset only
        #only select the subset of genes from the learning data
        domain = orange.Domain([ldata.domain.attributes[ai] for ai in gs], ldata.domain.classVar)
        datao = orange.ExampleTable(domain, ldata)
        example = orange.Example(domain, ex) #domains need to be the same
      
        ocena = setSig_example_geneset(example, datao)
        pathways[name] = ocena
        
    return pathways

class SetSig(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        return setSig_example(self.learndata, example, self.genesets)

class SetSigLearner(object):

    def __call__(self, data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
        data, oknames, gsetsnum = selectGenesetsData(data, matcher, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        return SetSig(learndata=data, genesets=gsetsnum)

if __name__ == "__main__":

    data = Orange.data.Table("iris")
    gsets = obiGeneSets.collections({
        #"ALL": ['sepal length', 'sepal width', 'petal length', 'petal width'],
        "f3": ['sepal length', 'sepal width', 'petal length'],
        "l3": ['sepal width', 'petal length', 'petal width'],
        })

    fp = 120
    ldata = orange.ExampleTable(data.domain, data[:fp])
    tdata = orange.ExampleTable(data.domain, data[fp:])

    matcher = obiGene.matcher([])

    choosen_cv = ["Iris-setosa", "Iris-versicolor"]
    #ass = AssessLearner()(data, matcher, gsets, rankingf=AT_loessLearner())
    ass = AssessLearner()(data, matcher, gsets, minPart=0.0)
    #ass = MeanLearner()(data, matcher, gsets, default=False)
    #ass = CORGsLearner()(data, matcher, gsets)
    #ass = MedianLearner()(data, matcher, gsets)
    #ass = PLSLearner()(data, matcher, gsets, classValues=choosen_cv, minPart=0.0)
    #ass = SetSigLearner()(ldata, matcher, gsets, classValues=choosen_cv, minPart=0.0)
    #ass = PCALearner()(ldata, matcher, gsets, classValues=choosen_cv, minPart=0.0)
    #ass = GSALearner()(ldata, matcher, gsets, classValues=choosen_cv, minPart=0.0)

    ar = defaultdict(list)
    for d in (list(ldata) + list(tdata))[:5]:
        for a,b in ass(d).items():
            ar[a].append(b)

    def pp1(ar):
        ol =  sorted(ar.items())
        print '\n'.join([ a.id + ": " +str(b) for a,b in ol])

    pp1(ar)

