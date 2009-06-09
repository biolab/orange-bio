import stats, orange, numpy, statc

def mean(l):
    return float(sum(l))/len(l)

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

    Usege similar to MeasureAttribute*.

    Standard deviation used for now returns minmally 0.2*|mi|, where mi=0 is adjusted to mi=1
    (as in gsea implementation).

    Can work only on data with two classes. If there are multiple class, then
    relevant class values can be specified on object initialization.
    By default the relevant classes are first and second class value
    from the domain.
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

        def stdev(l):
            return statc.std(l)

        def mean(l):
            return statc.mean(l)

        def stdevm(l):
            m = mean(l)
            std = stdev(l)
            #return minmally 0.2*|mi|, where mi=0 is adjusted to mi=1
            return max(std, 0.2*abs(1.0 if m == 0 else m))

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[-1].value == value and not ex[i].isSpecial() ]

        exa = avWCVal(self.a)
        exb = avWCVal(self.b)

        try:
            rval = (mean(exa)-mean(exb))/(stdevm(exa)+stdevm(exb))
            return rval
        except:
            #return some "middle" value -
            #TODO rather throw exception? 
            return 0

class MA_t_test(object):
    def __init__(self, a=None, b=None, prob=False):
        self.a = a
        self.b = b
        self.prob = prob
    def __call__(self, i, data):
        cv = data.domain.classVar
        #print data.domain

        #for faster computation. to save dragging many attributes along
        dom2 = orange.Domain([data.domain[i]], data.domain.classVar)
        data = orange.ExampleTable(dom2, data)
        i = 0

        if self.a == None: self.a = cv.values[0]
        if self.b == None: self.b = cv.values[1]

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[cv] == value and not ex[i].isSpecial() ]

        exa = avWCVal(self.a)
        exb = avWCVal(self.b)

        try:
            t, prob = stats.lttest_ind(exa, exb)
            return prob if self.prob else t
        except:
            return 1.0 if self.prob else 0.0

class MA_fold_change(object):
    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b
    def __call__(self, i, data):
        cv = data.domain.classVar
        #print data.domain

        #for faster computation. to save dragging many attributes along
        dom2 = orange.Domain([data.domain[i]], data.domain.classVar)
        data = orange.ExampleTable(dom2, data)
        i = 0

        if self.a == None: self.a = cv.values[0]
        if self.b == None: self.b = cv.values[1]

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[cv] == value and not ex[i].isSpecial() ]

        exa = avWCVal(self.a)
        exb = avWCVal(self.b)

        try:
            return mean(exa)/mean(exb)
        except:
            return 1

class MA_anova(object):
    def __init__(self, prob=False):
        self.prob = prob
    def __call__(self, i, data):
        cv = data.domain.classVar
        #print data.domain

        #for faster computation. to save dragging many attributes along
        dom2 = orange.Domain([data.domain[i]], data.domain.classVar)
        data = orange.ExampleTable(dom2, data)
        i = 0

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[cv] == value and not ex[i].isSpecial() ]

        data = [avWCVal(val) for val in cv.values]

        try:
            f, prob = stats.lF_oneway(*tuple(data))
            return prob if self.prob else f
        except:
            return 1.0 if self.prob else 0.0

import numpy as np
import numpy.ma as ma

class ExpressionSignificance_Test(object):
    def __new__(cls, data, useAttributeLabels, **kwargs):
        self = object.__new__(cls)
        if kwargs:
            self.__init__(data, useAttributeLabels)
            return self.__call__(**kwargs)
        else:
            return self
    
    def __init__(self, data, useAttributeLabels=False):
        self.data = data
        self.useAttributeLabels = useAttributeLabels
        self.attr_labels, self.data_classes = self._data_info(data)
        self.attributes = [attr for attr in self.data.domain.attributes if attr.varType in [orange.VarTypes.Continuous, orange.VarTypes.Discrete]]
        self.classes = np.array(self.attr_labels if useAttributeLabels else self.data_classes)
        self.keys = range(len(data)) if useAttributeLabels else self.attributes
        self.array, _, _ = data.toNumpyMA()
        if self.useAttributeLabels:
            self.array = ma.transpose(self.array)
#        self.dim = 1 if useAttributeLabels else 0  
        self.dim = 0
        
    def _data_info(self, data):
        return [set(attr.attributes.values()) for attr in data.domain.attributes], [ex.getclass() for ex in data] if data.domain.classVar else [None]*len(data)

    def _separete(self, source, value):
        return ([i for i, v in enumerate(source) if v==value],
                [i for i, v in enumerate(source) if v!=value])
        
    def test_indices(self, target):
        if self.useAttributeLabels:
            ind1 = [i for i, cl in enumerate(self.classes) if cl >= set([target])]
            ind2 = [i for i, cl in enumerate(self.classes) if not cl >= set([target])]
        else:
            ind1 = ma.nonzero(self.classes == target)[0]
            ind2 = ma.nonzero(self.classes != target)[0]
        return ind1, ind2
    
    def __call__(self, target):
        raise NotImplementedError()
    
class ExpressionSignificance_TTest(ExpressionSignificance_Test):
    def __call__(self, target):
        ind1, ind2 = self.test_indices(target)
        t, pval = attest_ind(self.array[ind1, :], self.array[ind2, :], dim=self.dim)
        return zip(self.keys,  zip(t, pval))
        
class ExpressionSignificance_FoldChange(ExpressionSignificance_Test):
    def __call__(self, target):
        ind1, ind2 = self.test_indices(target)
        a1, a2 = self.array[ind1, :], self.array[ind2, :]
        fold = ma.mean(a1, self.dim)/ma.mean(a2, self.dim)
        return zip(self.keys, fold)
    
def attest_ind(a, b, dim=None):
    if dim==None:
        dim = 0
        a, b = ma.ravel(a), ma.ravel(b)
    x1, x2 = ma.mean(a, dim), ma.mean(b, dim)
    v1, v2 = ma.var(a, dim), ma.var(b, dim)
    n1, n2 = a.shape[dim], b.shape[dim]
    df = float(n1+n2-2)
    svar = ((n1-1)*v1+(n2-1)*v2) / df
    t = (x1-x2)/ma.sqrt(svar*(1.0/n1 + 1.0/n2))
#    probs = abetai(0.5*df,0.5,float(df)/(df+t*t))
    try:
#        prob = [statc.betai(0.5*df,0.5,float(df)/(df+tsq)) for tsq in list(t*t)]
        prob = [statc.betai(0.5*df,0.5,df/(df+tsq)) if tsq is not ma.masked and df/(df+tsq) <= 1.0 else ma.masked  for tsq in t*t]
    except :
        print "tsq:", tsq, type(tsq), (t*t).shape
        raise
    return t, prob