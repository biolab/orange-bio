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

        def stdev(l):
            return stats.stdev(l)

        def stdevm(l):
            m = mean(l)
            std = stdev(l)
            #return minmally 0.2*|mi|, where mi=0 is adjusted to mi=1
            return max(std, 0.2*abs(1.0 if m == 0 else m))

        def avWCVal(value):
            return [ex[i].value for ex in data if ex[cv] == value and not ex[i].isSpecial() ]

        exa = avWCVal(self.a)
        exb = avWCVal(self.b)

        try:
            t, prob = stats.lttest_ind(exa, exb)
            return prob if self.prob else t
        except:
            return 1.0 if sefl.prob else 0.0

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
