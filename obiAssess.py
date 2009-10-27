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
import obiGene
from collections import defaultdict

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
       XR = XR - dot(numpy.array([t]).T, numpy.array([P[:,i]]))
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
    YStd = numpy.std(Y, axis = 0)
    XMean = numpy.mean(X, axis = 0)
    XStd = numpy.std(X, axis = 0)
    
    #FIXME: standard deviation should never be 0. Ask Lan, if the following
    #fix is ok.
    XStd = numpy.maximum(XStd, 10e-16)
    YStd = numpy.maximum(YStd, 10e-16)

    #X = (X-XMean)/XStd
    #Y = (Y-YMean)/YStd
    X = (X-XMean)
    Y = (Y-YMean)

    P = numpy.empty((mx,Ncomp))
    C = numpy.empty((my,Ncomp))
    T = numpy.empty((n,Ncomp))
    U = numpy.empty((n,Ncomp))
    B = numpy.zeros((Ncomp,Ncomp))
    W = numpy.empty((mx,Ncomp))
    E,F = X,Y

    dot = numpy.dot
    norm = numpy.linalg.norm

    #PLS1 - from Gutkin, shamir, Dror: SlimPLS

    for i in range(Ncomp):

        c = reduce(dot, [ F.T, E, E.T, F]) ** -0.5 #normalization factor
        w = c*dot(E.T,F) #w_a
        t = dot(E, w) #t_i
        p = dot(E.T, t)/dot(t.T, t) #p_i
        q = dot(F.T, t)/dot(t.T, t) #q_i
            
        E = E - dot(t, p.T)
        F = F - dot(t, q.T)

        T[:,i] = t.T
        W[:,i] = w.T
        P[:,i] = p.T
        

    """
    print "T", T

    #regenerate T - see if they match
    
    TR = numpy.empty((n,Ncomp))
    XR = X

    for i in range(Ncomp):
       t = dot(XR, W[:,i].T)
       XR = XR - dot(numpy.array([t]).T, numpy.array([P[:,i]]))
       TR[:,i] = t

    print "TR", TR
    """

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
    
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None, components=1):
        """
        If more that 1 components are used, _LC_componetsNumber is appended to 
        the name of the gene set.
        """

        data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
    
        constructt = {}

        #build weight matrices for every gene set
        for name, inds in gsetsnum.items():
            dom2 = orange.Domain([ data.domain.attributes[i1] for i1 in inds ], data.domain.classVar)
            data_gs = orange.ExampleTable(dom2, data)
            xmean, W, P, T = PLSCall(data_gs, nc=components, y=[data_gs.domain.classVar], save_partial=True)
            constructt[name] = inds, xmean, W, P

            """
            print "TO", T[0]
            pt2 = pls_transform(data[0], constructt[name])
            print "TR", pt2
            """

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

class SetSigOLD(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, example):
        #return a dictionary geneset: value for every sample
        test,_ = genesetsAsAttributes(self.learndata, [ example ], self.genesets, self.minGenes, self.maxGenes, domain=self.domain)
        #print test[0]
        return dict(  (at.name, test[0][at].value) for at in test.domain.attributes )

class SetSigOLDLearner(object):
    """
    Just applies a function taking attribute values of an example and
    produces to all gene sets.    
    """
    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
        #data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
        #    minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        gm = obiGene.matcher([obiGene.GMKEGG("hsa")])
        gm.set_targets([g.name for g in data.domain.attributes])
        genesets = genesetsInData(geneSets, gm)
        learn,domain = genesetsAsAttributes(data, data, genesets, minSize, maxSize)
        #learn here was the original output of tranformed learning set - we get the same output
        #if we call the SetSigOLD with the same example
        """
        print "LEARN"
        print learn.domain
        for ex in learn:
            print ex
        """
        return SetSigOLD(learndata=data, domain=domain, genesets=genesets, minGenes=minSize, maxGenes=maxSize)

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

    def __call__(self, data, organism, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
        data, oknames, gsetsnum = selectGenesetsData(data, organism, geneSets, \
            minSize=minSize, maxSize=maxSize, minPart=minPart, classValues=classValues)
        return SetSig(learndata=data, genesets=gsetsnum)

if __name__ == "__main__":
    
    data = orange.ExampleTable("DLBCL.tab")

    """
    data = orange.ExampleTable("sterolTalkHepa.tab")
    data = impute_missing(data)
    choosen_cv = [ "LK935_48h", "Rif_12h"]
    ncl = orange.EnumVariable("cl", values=choosen_cv)
    ncl.getValueFrom = lambda ex,rw: orange.Value(ncl, ex[-1].value)
    ndom = orange.Domain(data.domain.attributes, ncl)
    data = orange.ExampleTable(ndom, [ ex for ex in data if ex[-1].value in choosen_cv ])
    """

    choosen_cv = list(data.domain.classVar.values)

    fp = int(9*len(data)/10)

    ldata = orange.ExampleTable(data.domain, data[:fp])
    tdata = orange.ExampleTable(data.domain, data[fp:])


    gsets = obiGeneSets.collections(["steroltalk.gmt"], default=False)
    #gsets = obiGeneSets.collections(["C2.CP.gmt", "C5.MF.gmt", "C5.BP.gmt"], default=False)

    #ass = AssessLearner()(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False), rankingf=AT_loessLearner())
    #ass = MeanLearner()(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False))
    #ass = PLSLearner()(data, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False), classValues=choosen_cv)
    #ass = SetSigOLDLearner()(ldata, "hsa", obiGeneSets.collections(["steroltalk.gmt"], default=False), classValues=choosen_cv, minPart=0.0)
    ass = SetSigLearner()(ldata, "hsa", gsets, classValues=choosen_cv, minPart=0.0)

    ar = defaultdict(list)

    print data.domain.classVar.values

    for d in list(ldata) + list(tdata):
        for a,b in ass(d).items():
            ar[a].append(b)

    ol =  sorted(ar.items())
    print ol

    #print '\n'.join([ str(a) + ": " +str(b) for a,b in ol])
