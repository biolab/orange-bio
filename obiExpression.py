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
        return [set(attr.attributes.items()) for attr in data.domain.attributes], [ex.getclass() for ex in data] if data.domain.classVar else [None]*len(data)
        
    def test_indices(self, target, classes=None):
        classes = self.classes if classes is None else classes
        if self.useAttributeLabels:
            if type(target) in [list]:
                ind = [[i for i, cl in enumerate(self.classes) if cl >= set([t])] for t in target]
            else:
                ind1 = [i for i, cl in enumerate(self.classes) if cl >= set([target])]
                ind2 = [i for i, cl in enumerate(self.classes) if not cl >= set([target])]
                ind = [ind1, ind2]
        else:
            if type(target) in [list, tuple]:
                ind = [ma.nonzero(self.classes == t)[0] for t in target]
            else:
                ind1 = ma.nonzero(self.classes == target)[0]
                ind2 = ma.nonzero(self.classes != target)[0]
                ind = [ind1, ind2]
        return ind
    
    def __call__(self, target):
        raise NotImplementedError()
    
    def null_distribution(self, num, *args, **kwargs):
        kwargs = dict(kwargs)
        advance = lambda: None
        if "advance" in kwargs:
            advance = kwargs["advance"]
            del kwargs["advance"]
        results = []
        originalClasses = self.classes.copy()
        for i in range(num):
            np.random.shuffle(self.classes)
            results.append(self.__call__(*args, **kwargs))
            advance()
        self.classes = originalClasses
        return results
    
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
    
class ExpressionSignificance_SignalToNoise(ExpressionSignificance_Test):
    def __call__(self, target):
        ind1, ind2 = self.test_indices(target)
        a1, a2 = self.array[ind1, :], self.array[ind2, :]
        stn = (ma.mean(a1, self.dim) - ma.mean(a2, self.dim)) / (ma.sqrt(ma.var(a1, self.dim)) + ma.sqrt(ma.var(a2, self.dim)))
        return zip(self.keys, stn)
    
class ExpressionSignificance_ANOVA(ExpressionSignificance_Test):
    def __call__(self, target=None):
        if target is not None:
            indices = self.test_indices(target)
        else:
            indices = []
        f, prob = aF_oneway(*[self.array[ind, :] for ind in indices], **dict(dim=0))
        return zip(self.keys, zip(f, prob))
        
class ExpressionSignificance_ChiSquare(ExpressionSignificance_Test):
    def __call__(self, target):
        array = equi_n_discretization(self.array.copy(), intervals=5, dim=0)
        ind1, ind2 = self.test_indices(target)
        a1, a2 = array[ind1, :], array[ind2, :]
        dist1, dist2  = [], []
        dist = ma.zeros((array.shape[1], 2, 5))
        for i in range(5):
            dist1.append(ma.sum(ma.ones(a1.shape) * (a1 == i), 0))
            dist2.append(ma.sum(ma.ones(a2.shape) * (a2 == i), 0))
            dist[:, 0, i] = dist1[-1]
            dist[:, 1, i] = dist2[-1] 
        return zip(self.keys, achisquare_indtest(np.array(dist), dim=1))
        
class ExpressionSignificance_Info(ExpressionSignificance_Test):
    def __call__(self, target):
        array = equi_n_discretization(self.array.copy(), intervals=5, dim=1)
        
        ind1, ind2 = self.test_indices(target)
        a1, a2 = array[ind1, :], array[ind2, :]
        dist1, dist2 = [], []
        dist = ma.zeros((array.shape[1], 2, 5))
        for i in range(5):
            dist1.append(ma.sum(ma.ones(a1.shape) * (a1 == i), 0))
            dist2.append(ma.sum(ma.ones(a2.shape) * (a2 == i), 0))
            dist[:, 0, i] = dist1[-1]
            dist[:, 1, i] = dist2[-1]
        classinfo = entropy(np.array([len(ind1), len(ind2)]))
        E = ma.sum(entropy(dist, dim=1) * ma.sum(dist, 1), 1) / ma.sum(ma.sum(dist, 1), 1)
        return zip(self.keys, classinfo - E)
    
class ExpressionSignificance_MannWhitneyu(ExpressionSignificance_Test):
    def __call__(self, target):
        ind1, ind2 = self.test_indices(target)
        a, b = self.array[ind1, :], self.array[ind2, :]
#        results = [amannwhitneyu(a[:, i],b[:, i]) for i in range(a.shape[1])]
        results = [statc.mannwhitneyu(list(a[:, i]),list(b[:, i])) for i in range(a.shape[1])]
        
        return zip(self.keys, results)

def attest_ind(a, b, dim=None):
    """ Return the t-test statistics on arrays a and b over the dim axis.
    Returns both the t statistic as well as the p-value
    """
#    dim = a.ndim - 1 if dim is None else dim
    x1, x2 = ma.mean(a, dim), ma.mean(b, dim)
    v1, v2 = ma.var(a, dim), ma.var(b, dim)
    n1, n2 = (a.shape[dim], b.shape[dim]) if dim is not None else (a.size, b.size)
    df = float(n1+n2-2)
    svar = ((n1-1)*v1+(n2-1)*v2) / df
    t = (x1-x2)/ma.sqrt(svar*(1.0/n1 + 1.0/n2))
    if t.ndim == 0:
        return (t, statc.betai(0.5*df,0.5,df/(df+t**2)) if t is not ma.masked and df/(df+t**2) <= 1.0 else ma.masked)
    else:
        prob = [statc.betai(0.5*df,0.5,df/(df+tsq)) if tsq is not ma.masked and df/(df+tsq) <= 1.0 else ma.masked  for tsq in t*t]
        return t, prob

def aF_oneway(*args, **kwargs):
    dim = kwargs.get("dim", None)
    arrays = args
    means = [ma.mean(a, dim) for a in arrays]
    vars = [ma.var(a, dim) for a in arrays]
    lens = [ma.sum(ma.array(ma.ones(a.shape), mask=ma.asarray(a).mask), dim) for a in arrays]
    alldata = ma.concatenate(arrays, dim if dim is not None else 0)
    bign =  ma.sum(ma.array(ma.ones(alldata.shape), mask=alldata.mask), dim)
    sstot = ma.sum(alldata ** 2, dim) - (ma.sum(alldata, dim) ** 2) / bign
    ssbn = ma.sum([(ma.sum(a, dim) ** 2) / L for a, L in zip(arrays, lens)], dim)
#    print ma.sum(alldata, dim) ** 2 / bign, ssbn
    ssbn -= ma.sum(alldata, dim) ** 2 / bign
    sswn = sstot - ssbn
    dfbn = dfnum = float(len(args) - 1.0)
    dfwn = bign - len(args) # + 1.0
    F = (ssbn / dfbn) / (sswn / dfwn)
    if F.ndim == 0 and dfwn.ndim == 0:
        return (F,statc.betai(0.5 * dfwn, 0.5 * dfnum, dfwn/float(dfwn+dfnum*F)) if F is not ma.masked and dfwn/float(dfwn+dfnum*F) <= 1.0 \
                and dfwn/float(dfwn+dfnum*F) >= 0.0 else ma.masked)
    else:
        prob = [statc.betai(0.5 * dfden, 0.5 * dfnum, dfden/float(dfden+dfnum*f)) if f is not ma.masked and dfden/float(dfden+dfnum*f) <= 1.0 \
            and dfden/float(dfden+dfnum*f) >= 0.0 else ma.masked for dfden, f in zip (dfwn, F)]
        return F, prob
    
def achisquare_indtest(observed, dim=None):
    if observed.ndim == 2:
        observed = ma.array([observed])
        if dim is not None:
            dim += 1 
    if dim is None:
        dim = observed.ndim - 2
    rowtotal = ma.sum(observed, dim + 1)
    coltotal = ma.sum(observed, dim)
    total = ma.sum(rowtotal, dim)
    ones = ma.array(ma.ones(observed.shape))
    expected = ones * rowtotal.reshape(rowtotal.shape[:dim] + (-1, 1))
    a = ones * coltotal[..., np.zeros(observed.shape[dim], dtype=int),:]
    expected = expected * (a) / total.reshape((-1, 1, 1))
    chisq = ma.sum(ma.sum((observed - expected) ** 2 / expected, dim + 1), dim)
    return chisq
    
def equi_n_discretization(array, intervals=5, dim=1):
    count = ma.sum(ma.array(ma.ones(array.shape, dtype=int), mask=array.mask), dim)
    cut = ma.zeros(len(count), dtype=int)
    sarray = ma.sort(array, dim)
    r = count % intervals
    pointsshape = list(array.shape)
    pointsshape[dim] = 1
    points = []
    for i in range(intervals):
        cutend = cut + count / intervals + numpy.ones(len(r)) * (r > i)
        if dim == 1:
            p = sarray[range(len(cutend)), numpy.array(cutend, dtype=int) -1]
        else:
            p = sarray[numpy.array(cutend, dtype=int) -1, range(len(cutend))]
        points.append(p.reshape(pointsshape))
        cut = cutend
    darray = ma.array(ma.zeros(array.shape) - 1, mask=array.mask)
    darray[ma.nonzero(array <= points[0])] = 0
    for i in range(0, intervals):
        darray[ma.nonzero((array > points[i]))] = i + 1 
    return darray

def entropy(array, dim=None):
    if dim is None:
        array = array.ravel()
        dim = 0
    n = ma.sum(array, dim)
    array = ma.log(array) * array
    sum = ma.sum(array, dim)
    return (ma.log(n) - sum / n) / ma.log(2.0)
    