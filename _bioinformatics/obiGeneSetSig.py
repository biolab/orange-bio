from __future__ import absolute_import

import math
from collections import defaultdict

import scipy.stats

import numpy

import Orange, Orange.utils, statc

from .obiGsea import takeClasses
from .obiAssess import pca, PLSCall, corgs_activity_score
from . import obiExpression, obiGene, obiGeneSets, obiGsea, stats

class GeneSetTrans(object):

    __new__ = Orange.utils._orange__new__(object)

    def _mat_ni(self, data):
        """ With cached gene matchers. """
        if data.domain not in self._cache:
            self._cache[data.domain] = mat_ni(data, self.matcher)
        return self._cache[data.domain]

    def _match_instance(self, instance, geneset, takegenes=None):
        nm, name_ind = self._mat_ni(instance)
        genes = [ nm.umatch(gene) for gene in geneset ]
        if takegenes:
            genes = [ genes[i] for i in takegenes ]
        return nm, name_ind, genes

    def _match_data(self, data, geneset, odic=False):
        nm, name_ind = self._mat_ni(data)
        genes = [ nm.umatch(gene) for gene in geneset ]
        if odic:
            to_geneset = dict(zip(genes, geneset))
        takegenes = [ i for i,a in enumerate(genes) if a != None ]
        genes = [ genes[i] for i in takegenes ]
        if odic:
            return nm, name_ind, genes, takegenes, to_geneset
        else:
            return nm, name_ind, genes, takegenes

    def __init__(self, matcher=None, gene_sets=None, min_size=3, max_size=1000, min_part=0.1, class_values=None):
        self.matcher = matcher
        self.gene_sets = gene_sets
        self.min_size = min_size
        self.max_size = max_size
        self.min_part = min_part
        self.class_values = class_values
        self._cache = {}

    def __call__(self, data, weight_id=None):

        #selection of classes and gene sets
        data = takeClasses(data, classValues=self.class_values)
        nm,_ =  self._mat_ni(data)
        gene_sets = select_genesets(nm, self.gene_sets, self.min_size, self.max_size, self.min_part)

        #build a new domain
        newfeatures = self.build_features(data, gene_sets)
        newdomain = Orange.data.Domain(newfeatures, data.domain.class_var)
        return Orange.data.Table(newdomain, data)

    def build_features(self, data, gene_sets):
        return [ self.build_feature(data, gs) for gs in gene_sets ]


def normcdf(x, mi, st):
    return 0.5*(2. - stats.erfcc((x - mi)/(st*math.sqrt(2))))

class AT_edelmanParametric(object):

    def __init__(self, **kwargs):
        for a,b in kwargs.items():
            setattr(self, a, b)

    def __call__(self, nval):

        if self.mi1 == None or self.mi2 == None or self.st1 == None or self.st2 == None:
            return 0.0 

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
    Returns attribute transfromer for Edelman parametric measure for a
    given attribute in the dataset.  Edelman et al, 06. Modified a bit.
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

        val = nval

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
        try:
            ca = Orange.statistics.contingency.VarClass(data.domain.attributes[i], data)
            a =  Orange.statistics.estimate.ConditionalLoess(ca, nPoints=5)
            return AT_loess(condprob=a)
        except:
            return AT_loess(condprob=None)

def nth(l, n):
    return [a[n] for a in l]

class Assess(GeneSetTrans):
    """
    Uses the underlying GSEA code to select genes.
    Takes data and creates attribute transformations.
    """

    def __init__(self, rankingf=None, **kwargs):
        self.rankingf = rankingf
        if self.rankingf == None:
            self.rankingf = AT_edelmanParametricLearner()
        super(Assess, self).__init__(**kwargs)

    def build_features(self, data, gene_sets):

        attributes = []

        #attrans: { i_orig: ranking_function }
        attrans = [ self.rankingf(iat, data) for iat, at in enumerate(data.domain.attributes) ]

        nm_all, _ =  self._mat_ni(data)

        for gs in gene_sets:

            at = Orange.feature.Continuous(name=str(gs))

            geneset = list(gs.genes)
            nm, name_ind, genes, takegenes, to_geneset = self._match_data(data, geneset, odic=True)
            genes = set(genes)
            
            def t(ex, w, geneset=geneset, takegenes=takegenes, nm=nm, attrans=attrans):

                nm2, name_ind2, genes2 = self._match_instance(ex, geneset, takegenes)

                ex_atts = [ at.name for at in ex.domain.attributes ]
                new_atts = [ name_ind[nm.umatch(an)] if nm.umatch(an) != None else None
                    for an in ex_atts ]
                #new_atts: indices of genes in original data for that sample 
                #POSSIBLE REVERSE IMPLEMENTATION (slightly different
                #for data from different chips):
                #save pairs together and sort (or equiv. dictionary transformation)

                indexes = filter(lambda x: x[0] != None, zip(new_atts, range(len(ex_atts))))

                lcor = [ attrans[index_in_data](ex[index_in_ex].value) 
                    for index_in_data, index_in_ex in indexes if
                    ex[index_in_ex].value != '?' ]
                #indexes in original lcor, sorted from higher to lower values
                ordered = obiGsea.orderedPointersCorr(lcor) 
                #subset = list of indices, lcor = correlations, ordered = order
                subset = [ name_ind2[g] for g in genes2 ]
                return obiGsea.enrichmentScoreRanked(subset, lcor, ordered)[0] 

            at.get_value_from = t
            attributes.append(at)

        return attributes
   
def setSig_example_geneset(ex, data):
    """ Gets learning data and example with the same domain, both
    containing only genes from the gene set. """

    distances = [ [], [] ]    

    def pearsonr(v1, v2):
        return numpy.corrcoef([v1, v2])[0,1]

    def pearson(ex1, ex2):
        #leaves undefined elements out

        attrs = range(len(ex1.domain.attributes))
        vals1 = [ ex1[i].value for i in attrs ]
        vals2 = [ ex2[i].value for i in attrs ]

        common = [ True if v1 != "?" and v2 != "?" else False \
            for v1,v2 in zip(vals1,vals2) ]
        vals1 = [ v for v,c in zip(vals1, common) if c ]
        vals2 = [ v for v,c in zip(vals2, common) if c ]

        return numpy.corrcoef([vals1, vals2])[0,1]

    def ttest(ex1, ex2):
        try:
            return stats.lttest_ind(ex1, ex2)[0]
        except:
            return 0.0
    
    #maps class value to its index
    classValueMap = dict( [ (val,i) for i,val in enumerate(data.domain.class_var.values) ])
 
    #create distances to all learning data - save or other class
    for c in data:
        distances[classValueMap[c[-1].value]].append(pearson(c, ex))

    return ttest(distances[0], distances[1])

def mat_ni(data, matcher):
    nm = matcher([at.name for at in data.domain.attributes])
    name_ind = dict((n.name,i) for i,n in enumerate(data.domain.attributes))
    return nm, name_ind

def select_genesets(nm, gene_sets, min_size=3, max_size=1000, min_part=0.1):
    """ Returns a list of gene sets that have sizes in limits """

    def ok_sizes(gs):
        """compares sizes of genesets to limitations"""
        transl = filter(lambda x: x != None, [ nm.umatch(gene) for gene in gs.genes ])
        if len(transl) >= min_size \
            and len(transl) <= max_size \
            and float(len(transl))/len(gs.genes) >= min_part:
            return True
        return False

    return filter(ok_sizes, gene_sets) 

def vou(ex, gn, indices):
    """ returns the value or "?" for the given gene name gn"""
    if gn not in indices:
        return "?"
    else:
        return ex[indices[gn]].value

class SetSig(GeneSetTrans):

    def build_feature(self, data, gs):

        at = Orange.feature.Continuous(name=str(gs))

        def t(ex, w, gs=gs, data=data): #copy od the data
            geneset = list(gs.genes)

            nm, name_ind, genes, takegenes = self._match_data(data, geneset)
            nm2, name_ind2, genes2 = self._match_instance(ex, geneset, takegenes)

            domain = Orange.data.Domain([data.domain.attributes[name_ind[gene]] for gene in genes], data.domain.class_var)
            datao = Orange.data.Table(domain, data)
           
            #convert the example to the same domain
            exvalues = [ vou(ex, gn, name_ind2) for gn in genes2 ] + [ "?" ]
            example = Orange.data.Instance(domain, exvalues)

            return setSig_example_geneset(example, datao) #only this one is setsig specific
     
        at.get_value_from = t
        return at

class ParametrizedTransformation(GeneSetTrans):

    def _get_par(self, datao):
        pass
        
    def _use_par(self, ex, constructt):
        pass
    
    def build_feature(self, data, gs):

        at = Orange.feature.Continuous(name=str(gs))

        geneset = list(gs.genes)
        nm, name_ind, genes, takegenes = self._match_data(data, geneset)
        domain = Orange.data.Domain([data.domain.attributes[name_ind[gene]] for gene in genes], data.domain.class_var)
        datao = Orange.data.Table(domain, data)

        constructt = self._get_par(datao)

        def t(ex, w, geneset=geneset, constructt=constructt, takegenes=takegenes, domain=domain):
            nm2, name_ind2, genes2 = self._match_instance(ex, geneset, takegenes)
          
            #convert the example to the same domain
            exvalues = [ vou(ex, gn, name_ind2) for gn in genes2 ] + [ "?" ]
            ex = numpy.array(exvalues[:-1])

            return self._use_par(ex, constructt)
            
        at.get_value_from = t
        return at

class PLS(ParametrizedTransformation):

    def _get_par(self, datao):
        return PLSCall(datao, nc=1, y=[datao.domain.class_var])
        
    def _use_par(self, ex, constructt):
        xmean, W, P, _ = constructt
        ex = ex - xmean # same input transformation

        nc = W.shape[1]

        TR = numpy.empty((1, nc))
        XR = ex

        dot = numpy.dot

        for i in range(nc):
           t = dot(XR, W[:,i].T)
           XR = XR - t*numpy.array([P[:,i]])
           TR[:,i] = t

        return TR[0][0]
        
class PCA(ParametrizedTransformation):

    def _get_par(self, datao):
        return pca(datao)

    def _use_par(self, arr, constructt):
        evals, evect, xmean = constructt

        arr = arr - xmean # same input transformation - a row in a matrix
        ev0 = evect[0] #this is a row in a matrix - do a dot product
        a = numpy.dot(arr, ev0)

        return a

class SimpleFun(GeneSetTrans):

    def build_feature(self, data, gs):

        at = Orange.feature.Continuous(name=str(gs))

        def t(ex, w, gs=gs):
            geneset = list(gs.genes)
            nm2, name_ind2, genes2 = self._match_instance(ex, geneset)
           
            exvalues = [ vou(ex, gn, name_ind2) for gn in genes2 ] + [ "?" ]
            exvalues = filter(lambda x: x != "?", exvalues)

            return self.fn(exvalues)
     
        at.get_value_from = t
        return at

class Mean(SimpleFun):

    def __init__(self, **kwargs):
       self.fn = numpy.mean
       super(Mean, self).__init__(**kwargs)

class Median(SimpleFun):

    def __init__(self, **kwargs):
       self.fn = numpy.median
       super(Median, self).__init__(**kwargs)

class GSA(GeneSetTrans):

    def build_features(self, data, gene_sets):

        attributes = []

        def tscorec(data, at, cache=None):
            ma = obiExpression.MA_t_test()(at,data)
            return ma

        tscores = [ tscorec(data, at) for at in data.domain.attributes ]

        def to_z_score(t):
            return float(scipy.stats.norm.ppf(scipy.stats.t.cdf(t, len(data)-2)))

        zscores = map(to_z_score, tscores)

        for gs in gene_sets:

            at = Orange.feature.Continuous(name=str(gs))

            geneset = list(gs.genes)
            nm, name_ind, genes, takegenes, to_geneset = self._match_data(data, geneset, odic=True)
            #take each gene only once
            genes = set(genes)

            D = numpy.mean([max(zscores[name_ind[g]],0) for g in genes]) \
                + numpy.mean([min(zscores[name_ind[g]],0) for g in genes])

            if D >= 0:
                consider_genes = [ to_geneset[g] for g in genes if zscores[name_ind[g]] > 0.0 ]
            else:
                consider_genes = [ to_geneset[g] for g in genes if zscores[name_ind[g]] < 0.0 ]

            def t(ex, w, consider_genes=consider_genes):
                nm2, name_ind2, genes2 = self._match_instance(ex, consider_genes)
              
                #convert the example to the same domain
                exvalues = [ vou(ex, gn, name_ind2) for gn in genes2 ] + [ "?" ]
                exvalues = filter(lambda x: x != "?", exvalues)
              
                return numpy.mean(exvalues)

            at.get_value_from = t
            attributes.append(at)

        return attributes

def tscorec(data, at, cache=None):
    """ Cached attribute  tscore calculation """
    if cache != None and at in cache: return cache[at]
    ma = obiExpression.MA_t_test()(at,data)
    if cache != None: cache[at] = ma
    return ma

def nth(l, n):
    return [a[n] for a in l]

def compute_corg(data, inds, tscorecache):
    """
    Compute CORG for this geneset specified with gene inds
    in the example table. Output is the list of gene inds
    in CORG.

    """
    #order member genes by their t-scores: decreasing, if av(t-score) >= 0,
    #else increasing
    tscores = [ tscorec(data, at, tscorecache) for at in inds ]
    sortedinds = nth(sorted(zip(inds,tscores), key=lambda x: x[1], \
        reverse=numpy.mean(tscores) >= 0), 0)

    def S(corg):
        """ Activity score separation - S(G) in 
        the article """
        asv = Orange.feature.Continuous(name='AS')
        asv.getValueFrom = lambda ex,rw: Orange.data.Value(asv, corgs_activity_score(ex, corg))
        data2 = Orange.data.Table(Orange.data.Domain([asv], data.domain.classVar), data)
        return abs(tscorec(data2, 0)) #FIXME absolute - nothing in the article abs()
            
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

class CORGs(ParametrizedTransformation):
    """
    WARNING: input has to be z_ij table! each gene needs to be normalized
    (mean=0, stdev=1) for all samples.
    """

    def __call__(self, *args, **kwargs):
        self.tscorecache = {} #reset a cache
        return super(CORGs, self).__call__(*args, **kwargs)

    def build_feature(self, data, gs):

        at = Orange.feature.Continuous(name=str(gs))
        geneset = list(gs.genes)

        nm, name_ind, genes, takegenes, to_geneset = self._match_data(data, geneset, odic=True)
        indices = compute_corg(data, [ name_ind[g] for g in genes ], self.tscorecache)

        ind_names = dict( (a,b) for b,a in name_ind.items() )
        selected_genes = sorted(set([to_geneset[ind_names[i]] for i in indices]))
            
        def t(ex, w, corg=selected_genes): #copy od the data
            nm2, name_ind2, genes2 = self._match_instance(ex, corg, None)
            exvalues = [ vou(ex, gn, name_ind2) for gn in genes2 ]
            return sum(v if v != '?' else 0.0 for v in exvalues)/len(corg)**0.5
     
        at.get_value_from = t
        return at

if __name__ == "__main__":

    data = Orange.data.Table("iris")
    gsets = obiGeneSets.collections({
        #"ALL": ['sepal length', 'sepal width', 'petal length', 'petal width'],
        "f3": ['sepal length', 'sepal width', 'petal length'],
        "l3": ['sepal width', 'petal length', 'petal width'],
        })
    matcher = obiGene.matcher([])
    choosen_cv = ["Iris-setosa", "Iris-versicolor"]

    """
    data = Orange.data.Table("DLBCL_200a")
    gsets = obiGeneSets.collections((("KEGG",),"9606"))
    matcher = obiGene.matcher([obiGene.GMKEGG("hsa")])
    choosen_cv = None
    """

    def to_old_dic(d, data):
        ar = defaultdict(list)
        for ex1 in data:
            ex = d(ex1)
            for a,x in zip(d.attributes, ex):
                ar[a.name].append(x.value)
        return ar

    def pp2(ar):
        ol =  sorted(ar.items())
        print '\n'.join([ a + ": " +str(b) for a,b in ol])

    ass = Assess(data, matcher=matcher, gene_sets=gsets, class_values=choosen_cv, min_part=0.0)
    ar = to_old_dic(ass.domain, data[:5])
    pp2(ar)
