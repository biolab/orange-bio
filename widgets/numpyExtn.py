###################################################################################
## numpy extension by PJ (peter.juvan@fri.uni-lj.si)
##
## TODO: fix default axis asignment (in several funcions 0 is the dedault):
##      - for numpy default is None
##      - for numpy.me default is still -1 (will this be fixed in future releases?)
###################################################################################

import numpy

###################################################################################
## pretty print
###################################################################################

def strSimMtrx(m, names, numPos=None, sort=False):
    """ print similarity matrix with fixed number of positions,
    appends names, optionally sorts names and m
    """
    assert len(names) == m.shape[0]
    assert len(m.shape) == 2
    m = numpy.asarray(m)
    if sort == True:
        nameIndList = zip(names, range(len(names)))
        nameIndList.sort()
        names = map(lambda x: x[0], nameIndList)
        sortArgs = map(lambda x: x[1], nameIndList)
        m = m[sortArgs][:,sortArgs]
    out = ""
    if numPos == None:
        numPos = max(map(lambda n: len(n), names))
    # number of decimal places
    numWhole = len(str(int(numpy.max(m))))
    numDec = numPos - numWhole - 1
    if numDec < 0:
        numDec = 0
    # print header row
    printStr = "%%s %%%is" % numPos
    out += reduce(lambda a,b: printStr % (a, b[:numPos]), names, " "*numPos)
    # print other rows
    printStrName = "%%%is" % numPos
    printStrMtrx2 = "%%.%if" % numDec
    for rowIdx,row in enumerate(m):
        rowStr = map(lambda x: printStrMtrx2 % x, row)
        out += "\n%s" % reduce(lambda a,b: printStr % (a, b[:numPos]), rowStr, printStrName % names[rowIdx][:numPos])
    return out


###################################################################################
## spearman r / pearson r between values which are present in both arrays
## uses scipy.stats
###################################################################################

def spearmanrMA(m1, m2):
    # calculate spearmanr between values which are present in both arrays
    try:
        import scipy.stats
    except ImportError:
        print "scipy.stats not found; you can download it from www.scipy.org"
    nomask12 = numpy.logical_not(numpy.logical_or(numpy.ma.getmaskarray(m1), numpy.ma.getmaskarray(m2)))
    return scipy.stats.spearmanr(m1.compress(nomask12), m2.compress(nomask12))

def pearsonrMA(m1, m2):
    # calculate pearson between values which are present in both arrays
    try:
        import scipy.stats
    except ImportError:
        print "scipy.stats not found; you can download it from www.scipy.org"
    nomask12 = numpy.logical_not(numpy.logical_or(numpy.ma.getmaskarray(m1), numpy.ma.getmaskarray(m2)))
    return scipy.stats.pearsonr(m1.compress(nomask12), m2.compress(nomask12))


###################################################################################
## argmaxMA(m,axis=-1)
## replacement for numpy.ma.argmax which does not work with masked arrays:
##  numpy.ma.argmax(numpy.ma.asarray([1,2])*numpy.ma.masked) => 0
##  numpyextn.argmax(numpy.ma.asarray([1,2])*numpy.ma.masked) => numpy.ma.masked
###################################################################################

def argmaxMA(m,axis=-1):
    return numpy.ma.where(numpy.all(numpy.ma.getmaskarray(m),axis), numpy.ma.masked, numpy.ma.argmax(m,axis))


###################################################################################
## write to tab delimited file, replace masked data by "?" by default
###################################################################################

def tofileTabMa(m, fileName, fill="?"):
    assert numpy.rank(m) == 2
    m = numpy.ma.asarray(m)
    ms = numpy.asarray(numpy.ma.filled(m), "str")
    mss = numpy.where(numpy.ma.getmaskarray(m), fill, ms)
    fd = file(fileName, "w")
    for row in mss:
        fd.write("%s\n" % reduce(lambda a,b: "%s\t%s" % (a,b), row))
    fd.close()


###################################################################################
## median
###################################################################################

def median(m,axis=0):
    """Returns median of the given nonmasked m array along the given axis.
    """
    if numpy.ma.isMA(m):
       raise TypeError, "wrong type for m (%s), try using medianMA instead" % type(m).__name__
    m = numpy.asarray(m, numpy.float)
    # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    m = numpy.transpose(m, [axis]+range(axis)+range(axis+1,numpy.rank(m)))
    return numpy.median(m)


def medianMA(m,axis=0):
    """Returns median of the given masked array m along the given axis.
    """
    return percentilesMa(m,0.5,axis)


def percentilesMa(m,perc,axis=0):
    """Returns the percentiles of the given masked array m along the given axis.
    """
    assert 0 < perc < 1
    m = numpy.ma.asarray(m, numpy.float)
    # 2008-06-23: BUG: does not work with negative axis values: 0 axis is always taken
    # assert -len(m.shape) <= axis < len(m.shape)
    assert 0 <= axis < len(m.shape)
    # 2008-06-23 BUG: at arbitrary axis, a single element is just returned
    ##    if len(m.shape) == 0 or (len(m.shape)==1 and m.shape[0]==1):
    if m.shape[axis] == 1:
        return m[axis]
    elif len(m.shape) == 1:
        # returns (1,) array or ma.masked
        mCount = numpy.ma.count(m)
        if mCount == 0:
            return numpy.ma.masked
        elif mCount == 1:
            return m.compressed()[0]
        else:
            _k = float(perc)*(mCount-1)
            k = int(_k)
            d = _k - k
            mIndSort = numpy.ma.argsort(m,fill_value=1e20)
            return m[mIndSort[k]] + d*(m[mIndSort[k+1]]-m[mIndSort[k]])
    else:
        # count the number of nonmasked indices along the given axis
        _k = float(perc) * (numpy.ma.count(m, axis=axis)-1)
        k = _k.astype(numpy.int)
        d = _k - k
        # prepare indices for other axis (except for the given)
        cind = numpy.ma.indices(numpy.ma.shape(k))
        # indices of m sorted by increasing value at axis 0; masked values placed at the end;
        mtIndSort = numpy.ma.argsort(m,axis,fill_value=1e20)
        # get indices to address median elements from m
        # tuple(takeInd1): all lists in the tuple must be of the same shape,
        #   the result mtIndSort[tuple(takeInd1)] is also of that shape
        takeInd1 = cind.tolist()
        takeInd1.insert(axis,k.tolist())
        medInd1 = mtIndSort[tuple(takeInd1)]
        takeInd2 = cind.tolist()
        takeInd2.insert(axis,(k+1).tolist())
        medInd2 = mtIndSort[tuple(takeInd2)]
        # get both median elements; if c[i] odd: med1[i] == med2[i]
        takeMed1 = cind.tolist()
        takeMed1.insert(axis,medInd1.tolist())
        med1 = m[tuple(takeMed1)]
        takeMed2 = cind.tolist()
        takeMed2.insert(axis,medInd2.tolist())
        med2 = m[tuple(takeMed2)]
##        if __name__=="__main__":
##            print "m\n",numpy.ma.filled(m,-1)
##            print "k", k
##            print "[(k-1)/2,k/2]", [(k-1)/2,k/2]
##            print "cind\n",cind
##            print "mtIndSort\n",mtIndSort
##            print "medInd1\n",medInd1
##            print "medInd2\n",medInd2
##            print "med1\n",med1
##            print "med2\n",med2
        return med1 + numpy.ma.filled(d*(med2-med1),0)

if __name__ == "__main__":
    # 2008-06-23: bug with negative indices (not fixed, added assertion
##    print "first +", medianMA(numpy.asarray([[10,1],[2,11]]),0)
##    print "first -", medianMA(numpy.asarray([[10,1],[2,11]]),-2)
##    print "second +", medianMA(numpy.asarray([[10,1],[2,11]]),1)
##    print "second -", medianMA(numpy.asarray([[10,1],[2,11]]),-1)
    # 2008-06-23: bug with medianMA(numpy.asarray([[149.0, 87.0]]), 0)
    print medianMA(numpy.asarray([149.0, 87.0]), 0)
    print medianMA(numpy.asarray([[149.0, 87.0]]), 0)
    print medianMA(numpy.asarray([[149.0, 87.0]]), 1)


## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by PJ

import math
import numpy.oldnumeric as Numeric
import numpy.oldnumeric.ma as MA
import numpy.oldnumeric.mlab as MLab, numpy.oldnumeric.linear_algebra as LinearAlgebra


#####################################################################################
## LOWESS (from Biopython, extended to interpolation/extrapolation)
#####################################################################################

def lowess2(x, y, xest, f=2./3., iter=3):
    """Returns estimated values of y in data points xest (or None if estimation fails).
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations."""
    x = Numeric.asarray(x, 'd')
    y = Numeric.asarray(y, 'd')
    xest = Numeric.asarray(xest, 'd')
    n = len(x)
    nest = len(xest)
    r = min(int(Numeric.ceil(f*n)),n-1) # radius: num. of points to take into LR
    h = [Numeric.sort(abs(x-x[i]))[r] for i in range(n)]    # distance of the r-th point from x[i]
    w = Numeric.clip(abs(([x]-Numeric.transpose([x]))/h),0.0,1.0)
    w = 1-w*w*w
    w = w*w*w
    hest = [Numeric.sort(abs(x-xest[i]))[r] for i in range(nest)]    # r-th min. distance from xest[i] to x
    west = Numeric.clip(abs(([xest]-Numeric.transpose([x]))/hest),0.0,1.0)  # shape: (len(x), len(xest)
    west = 1-west*west*west
    west = west*west*west
    yest = Numeric.zeros(n,'d')
    yest2 = Numeric.zeros(nest,'d')
    delta = Numeric.ones(n,'d')
    try:
        for iteration in range(iter):
            # fit xest
            for i in range(nest):
                weights = delta * west[:,i]
                b = Numeric.array([sum(weights*y), sum(weights*y*x)])
                A = Numeric.array([[sum(weights), sum(weights*x)], [sum(weights*x), sum(weights*x*x)]])
                beta = LinearAlgebra.solve_linear_equations(A,b)
                yest2[i] = beta[0] + beta[1]*xest[i]
            # fit x (to calculate residuals and delta)
            for i in range(n):
                weights = delta * w[:,i]
                b = Numeric.array([sum(weights*y), sum(weights*y*x)])
                A = Numeric.array([[sum(weights), sum(weights*x)], [sum(weights*x), sum(weights*x*x)]])
                beta = LinearAlgebra.solve_linear_equations(A,b)
                yest[i] = beta[0] + beta[1]*x[i]
            residuals = y-yest
            s = MLab.median(abs(residuals))
            delta = Numeric.clip(residuals/(6*s),-1,1)
            delta = 1-delta*delta
            delta = delta*delta
    except LinearAlgebra.LinAlgError:
        print "Warning: NumExtn.lowess2: LinearAlgebra.solve_linear_equations: Singular matrix"
        yest2 = None
    return yest2

##x = Numeric.array([0,2.1,4.5, 6.], 'd')
##xs = Numeric.array([0,1,2,3,4,5,6], 'd')
##y = Numeric.array([0.8, 2.3, 2.9, 3.1], 'd')
##f = 2./3.
##iter=3
##import Bio.Statistics.lowess
##print Bio.Statistics.lowess.lowess(x, y, f, iter)
##print lowess2(x, y, xs, f, iter)

def lowessW(x, y, xest, f=2./3., iter=3, dWeights=None, callback=None):
    """Returns estimated values of y in data points xest (or None if estimation fails).
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.

    Data points may be assigned weights; if None, all weights equal 1.
    """
    x = Numeric.asarray(x, 'd')
    y = Numeric.asarray(y, 'd')
    xest = Numeric.asarray(xest, 'd')
    n = len(x)
    if n <> len(y):
        raise AttributeError, "Error: lowessW(x,y,xest,f,iter,dWeights): len(x)=%i not equal to len(y)=%i" % (len(x), len(y))
    nest = len(xest)
    # weights of data points (optional)
    if dWeights <> None:
        dWeights = Numeric.asarray(dWeights, 'd')
        if len(dWeights) <> n:
            raise AttributeError, "Error: lowessW(x,y,xest,f,iter,dWeights): len(dWeights)=%i not equal to len(x)=%i" % (len(dWeights), len(x))
##        dWeights = dWeights.reshape((n,1))
    else:
##        dWeights = Numeric.ones((n,1))
        dWeights = Numeric.ones((n,))
    r = min(int(Numeric.ceil(f*n)),n-1) # radius: num. of points to take into LR
    h = [Numeric.sort(abs(x-x[i]))[r] for i in range(n)]    # distance of the r-th point from x[i]
    w = Numeric.clip(abs(([x]-Numeric.transpose([x]))/h),0.0,1.0)
    w = 1-w*w*w
    w = w*w*w
    hest = [Numeric.sort(abs(x-xest[i]))[r] for i in range(nest)]    # r-th min. distance from xest[i] to x
    west = Numeric.clip(abs(([xest]-Numeric.transpose([x]))/hest),0.0,1.0)  # shape: (len(x), len(xest))
    west = 1-west*west*west
    west = west*west*west
    yest = Numeric.zeros(n,'d')
    yest2 = Numeric.zeros(nest,'d')
    delta = Numeric.ones(n,'d')
    try:
        for iteration in range(int(iter)):
            # fit xest
            for i in range(nest):
##                print delta.shape, west[:,i].shape, dWeights.shape
                weights = delta * west[:,i] * dWeights
                b = Numeric.array([sum(weights*y), sum(weights*y*x)])
                A = Numeric.array([[sum(weights), sum(weights*x)], [sum(weights*x), sum(weights*x*x)]])
                beta = LinearAlgebra.solve_linear_equations(A,b)
                yest2[i] = beta[0] + beta[1]*xest[i]
            # fit x (to calculate residuals and delta)
            for i in range(n):
                weights = delta * w[:,i] * dWeights
                b = Numeric.array([sum(weights*y), sum(weights*y*x)])
                A = Numeric.array([[sum(weights), sum(weights*x)], [sum(weights*x), sum(weights*x*x)]])
                beta = LinearAlgebra.solve_linear_equations(A,b)
                yest[i] = beta[0] + beta[1]*x[i]
            residuals = y-yest
            s = MLab.median(abs(residuals))
            delta = Numeric.clip(residuals/(6*s),-1,1)
            delta = 1-delta*delta
            delta = delta*delta
            if callback: callback()
    except LinearAlgebra.LinAlgError:
        print "Warning: NumExtn.lowessW: LinearAlgebra.solve_linear_equations: Singular matrix"
        yest2 = None
    return yest2


##if __name__ == "__main__":
##    x1 = numpy.asarray([0,1,2,3,4,5,6,7,8,9])
##    xe1 = numpy.asarray([-1,0,1,2,3,4,5,6,7,8,9,10])
##    xe2 = numpy.asarray([4,5,6])
##    y1 = numpy.asarray([0.1, 0.9, 2.05, 3.11, 3.99, 4.95, 5.88, 6.5, 6.8, 7])
##    w1 = numpy.asarray([1,   0.1,1,    1,    0.1, 0.1, 0.1, 0.1,0.1,0.1])
##
##    l2 = lowess2(x1, y1, xe1, f=0.3, iter=3)
##    print l2
##    l3 = lowessW(x1, y1, xe2, f=0.3, iter=3, dWeights=w1)
##    print l3 #, "\n", w3, "\n", west3



#####################################################################################
#####################################################################################

def denary2Binary(n):
    '''convert denary integer n to binary string bStr'''
    bStr = ''
    if n < 0:  raise ValueError, "must be a positive integer"
    if n == 0: return '0'
    while n > 0:
        bStr = str(n % 2) + bStr
        n = n >> 1
    return bStr


#####################################################################################
## position2index:
##  given the position of an element in a flat array
##  get the indices of that element in an array of given shape
##
## getIndices(m,val):
##  Input: m=[[1,2],[2,3]], val=2
##  Output: Numeric.array([(0,1),(1,0)])
#####################################################################################

def position2index(n, shape):
    """Given a position in a flat array (n), returns indices in an array of given shape.
    """
    raise DeprecationWarning, "Depricated, use positions2indices(pos, shape) instead."
    if n < 0:  raise ValueError, "must be a positive integer"
    shapeL = list(shape)
    shapeL.reverse()
    ind = []
    for shp in shapeL:
        ind = [n % shp] + ind
        n = int(n / shp)
    return ind


def positions2indices(pos, shape):
    """Given postions and shape, returns indices shaped (numPositions, len(shape));
    e.g.: pos=[0,2,8], shape=(3,3) ->  [[0,0],[2,0],[2,2]]
    """
    lenShape = len(shape)
    pos = numpy.asarray(pos)
    if pos.shape[0] == 0:
        return numpy.empty([0] + list(shape[1:]))
    posMax = numpy.max(pos)
    assert posMax < numpy.multiply.reduce(shape), "Error: Position too large for a given shape."
    assert numpy.min(pos) >= 0, "Error, position cannot be negative."
    ind = numpy.ones((pos.shape[0], lenShape))
    for dim in range(lenShape-1,-1,-1):
        ind[:,dim] = numpy.mod(pos, shape[dim])
        pos = numpy.divide(pos, shape[dim])
    return numpy.asarray(ind, numpy.int)


def indices2positions(ind, shape):
    """Given indices shaped (numPositions, len(shape)), returns corresponding positions in a falt array; output shape = (ind.shape[0],);
    e.g.: ind=[[0,0],[2,0],[2,2]] , shape=(3,3) -> [0,6,8]
    """
    ind = numpy.asarray(ind)
    if ind.shape[0] == 0:
        return numpy.empty((0,))
    assert len(ind.shape) == 2 and ind.shape[1] == len(shape), "Error: ind should be 2D array, ind.shape[1] should match len(shape)."
    assert numpy.less(numpy.max(ind,0),numpy.asarray(shape)).all(), "Error: indices do not fit the shape."
    pos = numpy.zeros((ind.shape[0],))
    for i,shp in enumerate(shape):
        pos += ind[:,i]*numpy.multiply.reduce(shape[i+1:])
    return numpy.asarray(pos,numpy.int)
    

def getPositions(m, val):
    """Input: arbitrary (masked) array and a value from that array;
    Output: array of positions of the given value in a flat m;
    """
    m = MA.asarray(m)
    return Numeric.compress(MA.equal(MA.ravel(m),val), Numeric.arange(Numeric.multiply.reduce(m.shape)))

def getIndices(m, val):
    """Input: arbitrary (masked) array and a value from that array;
    Output: array of indices corresponding to positions of the given value in array m;
    Output shape: (numb. of val values in m, len(m.shape)).
    """
    posFlat = getPositions(m, val)
    return positions2indices(posFlat, m.shape)


#####################################################################################
## Indices <-> Condition
#####################################################################################

def condition2indices(condition):
    """Input: condition=[1,0,0,1]; output: indices=[0,3]
    """
    condition = Numeric.asarray(condition)
    assert len(condition.shape) == 1
    return Numeric.compress(condition, Numeric.arange(condition.shape[0]))


def indices2condition(indices, length):
    """Input: indices=[0,3]; output: condition=[1,0,0,1]
    """
    indices = Numeric.asarray(indices)
    assert len(indices.shape) == 1
    assert length >= indices.shape[0]
    c = Numeric.zeros((length,), Numeric.Int)
    Numeric.put(c, indices, 1)
    return c

#####################################################################################
## Inverse permutation
##  permutInverse([0,2,3,5,4,1]) -> [0,5,1,2,4,3,]
##  permutInverse([0,5,1,2,4,3,]) -> [0,2,3,5,4,1]
#####################################################################################

#from Permutation import *

def permutInverse(n):
    """Returns inverse permutation given integers in range(len(n)),
    such that permitInverse(permutInverse(range(4)))==range(4).
    """
    n = Numeric.asarray(n)
    pInv = Numeric.argsort(n)
    assert Numeric.all(Numeric.equal(n, Numeric.argsort(pInv))), "Inverse not successful; input should be permutation of range(len(input))."
    return pInv

###################################################################################################################################
## Generator of random matrix of permutations
## where each row and each column represents a permutation of n elements.
###################################################################################################################################

def getRandomMtrxPermut(n):
    """Generator of random matrix of permutations;
    input: number of elements;
    returns a matrix where each row and each column represents a permutation of n elements.
    """
    mx = numpy.zeros((n,n), numpy.int)
    pArr = numpy.random.permutation(n)
    for i in range(n):
        mx[i] = numpy.concatenate((pArr[i:], pArr[:i]))
    return numpy.take(mx,numpy.random.permutation(n),0)   


#####################################################################################
## Ranks:
##      - Numeric (instead of statc.rankdata), in range 1 ... shape[0]
##      - MA: masked elements are not ranked; range: 1 ... #non-masked_values
#####################################################################################

def rankData(n, inverse=False):
    """Returns ranks of 1D Numeric array in range 1...shape[0].
    """
    n = Numeric.asarray(n)
    assert Numeric.rank(n) == 1
    r = Numeric.zeros(n.shape[0], Numeric.Float)
    Numeric.put(r, Numeric.argsort(n), Numeric.arange(n.shape[0]))
    if inverse:
        return -1*r+n.shape[0]
    else:
        return r+1

def rankDataMA(m, inverse=False):
    """Returns ranks of 1D masked array; masked values ignored, range 1...#non-masked_values.
    """
    m = MA.asarray(m)
    assert MA.rank(m) == 1
    fill_val = m.fill_value()
    m.set_fill_value(MA.maximum(m) + 1)
    r = MA.zeros(m.shape[0], Numeric.Float)
    MA.put(r, MA.argsort(m), Numeric.arange(m.shape[0]))
    m.set_fill_value(fill_val)
    r = MA.array(r, mask=MA.getmaskarray(m))
    if inverse:
        return -1*r+MA.count(m)
    else:
        return r+1
    


#####################################################################################
## Unary functions: a mask is set only in places where both a arrays are masked.
##  e.g.:
##    minus_unary(10,--)  -> 10
##    minus_unary(--,--) -> --
#####################################################################################

def subtract_unary(a, b):
    """Returns a-b with masked values only in places where both a and b are masked.
    """
    a = MA.asarray(a)
    b = MA.asarray(b)
    el = MA.subtract(a.filled(0), b.filled(0))
    mask = Numeric.logical_and(MA.getmaskarray(a), MA.getmaskarray(b))
    return MA.array(el, mask=mask)

def add_unary(a, b):
    """Returns a+b with masked values only in places where both a and b are masked.
    """
    a = MA.asarray(a)
    b = MA.asarray(b)
    el = MA.add(a.filled(0), b.filled(0))
    mask = Numeric.logical_and(MA.getmaskarray(a), MA.getmaskarray(b))
    return MA.array(el, mask=mask)

def multiply_unary(a, b):
    """Returns a*b with masked values only in places where both a and b are masked.
    """
    a = MA.asarray(a)
    b = MA.asarray(b)
    el = MA.multiply(a.filled(1), b.filled(1))
    mask = Numeric.logical_and(MA.getmaskarray(a), MA.getmaskarray(b))
    return MA.array(el, mask=mask)

def divide_unary(a, b):
    """Returns a*b with masked values only in places where both a and b are masked.
    """
    a = MA.asarray(a)
    b = MA.asarray(b)
    el = MA.divide(a.filled(1), b.filled(1))
    mask = Numeric.logical_and(MA.getmaskarray(a), MA.getmaskarray(b))
    return MA.array(el, mask=mask)


#####################################################################################
## Logical OR / AND of two masked arrays where:
##  - the unary OR equals to the element: OR(0) == 0, OR(1) == 1
##  - the mask equals to logical AND of the corresponding masks
##
## For example:
##  m1 = MA.array([0,0,0,-,-,-,1,1,1], mask=[0,0,0,1,1,1,0,0,0])
##  m2 = MA.array([0,-,1,0,-,1,0,-,1], mask=[0,1,0,0,1,0,0,1,0])
##
## logical_unary_or(m1,m2) =
##       MA.array([0,0,1,0,-,1,1,1,1], mask=[0,0,0,0,1,0,0,0,0])
## The built-in logical_or returns:
##       MA.array([0,-,0,-,-,-,0,-,1], mask=[0,1,0,1,1,1,0,1,0])
##
## logical_unary_and(m1,m2) =
##       MA.array([0,0,0,0,-,1,0,1,1], mask=[0,0,0,0,1,0,0,0,0])
## The built-in logical_and returns:
##       MA.array([0,-,0,-,-,-,0,-,1], mask=[0,1,0,1,1,1,0,1,0])
#####################################################################################

def logical_unary_or(m1, m2):
    el = Numeric.logical_or(m1.filled(0), m2.filled(0))
    mask = Numeric.logical_and(MA.getmaskarray(m1), MA.getmaskarray(m2))
    return MA.array(el, mask=mask)    


def logical_unary_and(m1, m2):
    el = Numeric.logical_and(m1.filled(1), m2.filled(1))
    mask = Numeric.logical_and(MA.getmaskarray(m1), MA.getmaskarray(m2))
    return MA.array(el, mask=mask)    


#####################################################################################
## MA.dot fix:
##  - MA.dot returns does not set the masked values (it places 0 instead)
##  - MA.dot(x,y)[i,j] == MA.masked iff MA.add.reduce(x[i,:]*y[:,j]) == MA.masked
##  - done faste enugh - it requires two matrix multiplications
#####################################################################################

def dotMA(a, b):
    """Returns dot-product for MA arrays; fixed masked values.
    """
    a = MA.asarray(a)
    b = MA.asarray(b)
    ab = MA.dot(a,b)
    # fix masked values in ab (MA.dot returns 0 instead of MA.masked)
    nonMasked = Numeric.dot(1-MA.getmaskarray(a).astype(Numeric.Int), 1-MA.getmaskarray(b).astype(Numeric.Int))
    return MA.where(Numeric.equal(nonMasked,0), MA.masked, ab)

#####################################################################################
## Apply a binary function on a sequence of arrays
#####################################################################################

def apply_binary(aFunction, *args):
    assert len(args) > 1, "at least two arguments required"
    return reduce(lambda x,y: aFunction(x,y), list(args[1:]), args[0])
##    return reduce(lambda x,y: aFunction(x,y), list(args))
        

#####################################################################################
## Strictly upper (lower) triangular matrix <-> 1D array
#####################################################################################

def triangularPut(m1d, upper=1, lower=0):
    """Returns 2D masked array with elements of the given 1D array in the strictly upper (lower) triangle.
    Elements of the 1D array should be ordered according to the upper triangular part of the 2D matrix.
    The lower triangular part (if requested) equals to the transposed upper triangular part.
    If upper == lower == 1 a symetric matrix is returned.
    """
    assert upper in [0,1] and lower in [0,1], "[0|1] expected for upper / lower"
    m1d = MA.asarray(m1d)
    assert MA.rank(m1d) == 1, "1D masked array expected"
    m2dShape0 = math.ceil(math.sqrt(2*m1d.shape[0]))
    assert m1d.shape[0] == m2dShape0*(m2dShape0-1)/2, "the length of m1d does not correspond to n(n-1)/2"
    if upper:
        if lower:
            mask = Numeric.fromfunction(lambda i,j: i==j, (m2dShape0, m2dShape0))
        else:
            mask = Numeric.fromfunction(lambda i,j: i>=j, (m2dShape0, m2dShape0))
    else:
        if lower:
            mask = Numeric.fromfunction(lambda i,j: i<=j, (m2dShape0, m2dShape0))
        else:
            mask = Numeric.ones((m2dShape0, m2dShape0))

    m2d = MA.ravel(MA.zeros((m2dShape0, m2dShape0), m1d.dtype.char))
    condUpperTriang = Numeric.fromfunction(lambda i,j: i<j, (m2dShape0, m2dShape0))
    putIndices = Numeric.compress(Numeric.ravel(condUpperTriang), Numeric.arange(0, m2dShape0**2, typecode=Numeric.Int))
    MA.put(m2d, putIndices, m1d)
    m2d = MA.reshape(m2d, (m2dShape0, m2dShape0))
    m2d = MA.where(condUpperTriang, m2d, MA.transpose(m2d))
    return MA.array(m2d, mask=Numeric.logical_or(mask, MA.getmaskarray(m2d)))


def triangularGet(m2d, upper=1):
    """Returns 1D masked array with elements from the upper (lower) triangular part of the given matrix.
    For a symetric matrix triangularGet(m2d, 0) and triangularGet(m2d, 1) return elements in different order.
    """
    assert upper in [0,1], "upper: [0|1] expected"
    m2d = MA.asarray(m2d)
    assert MA.rank(m2d) == 2, "2D (masked) array expected"
    if upper:
        takeInd = Numeric.compress(Numeric.ravel(Numeric.fromfunction(lambda i,j: i<j, m2d.shape)), Numeric.arange(0, Numeric.multiply.reduce(m2d.shape), typecode=Numeric.Int))
    else:
        takeInd = Numeric.compress(Numeric.ravel(Numeric.fromfunction(lambda i,j: i>j, m2d.shape)), Numeric.arange(0, Numeric.multiply.reduce(m2d.shape), typecode=Numeric.Int))
    return MA.ravel(m2d).take(takeInd)


#####################################################################################
## Put 1D array on the diagonal of the given 2D array and returns a copy.
#####################################################################################

def diagonalPut(m1d, m2d):
    """Puts the given 1D masked array into the diagonal of the given 2D masked array and returns a new copy of the 2D array.
    """
    m1d = MA.asarray(m1d)
    m2d = MA.asarray(m2d)
    assert MA.rank(m1d) == 1 and MA.rank(m2d) == 2, "1D and 2D masked array expected"
    assert m1d.shape[0] == m2d.shape[0] == m2d.shape[1], "the shape of the given arrays does not match"
    putIndices = Numeric.compress(Numeric.ravel(Numeric.fromfunction(lambda i,j: i==j, m2d.shape)), Numeric.arange(0, Numeric.multiply.reduce(m2d.shape), typecode=Numeric.Int))
    m2dShape = m2d.shape
    m2d = MA.ravel(m2d)
    MA.put(m2d, putIndices, m1d)
    return MA.reshape(m2d, m2dShape)
    

#####################################################################################
#### linear interpolation of missing values by the given axis
#####################################################################################
##
##def fillLinear1(ma, axis=0):
##    """BUG: does not handle successive missing values
##    """
##    maflat = MA.array(MA.ravel(ma))
##    mod = ma.shape[axis]
##    maskInd = Numeric.compress(MA.getmaskarray(maflat), Numeric.arange(maflat.shape[0]))
##    for idx in maskInd:
##        if (idx % mod) == 0:  # the first one
##            maflat[idx] = maflat[idx+1]
##        elif (idx % mod) == mod-1:
##            maflat[idx] = maflat[idx-1]
##        else:
##            maflat[idx] = (maflat[idx-1] + maflat[idx+1]) / 2.
##    return MA.reshape(maflat, ma.shape)
##
##
##def fillLinear2(ma, axis=0):
##    """BUG: does not handle successive missing values
##    """
##    maflat = MA.array(MA.ravel(ma))
##    mod = ma.shape[axis]
##    maskInd = Numeric.compress(MA.getmaskarray(maflat), Numeric.arange(maflat.shape[0]))
##    for idx in maskInd:
##        idxFirst = idx-(idx%mod)
##        idxNext = idxFirst+mod
##        if MA.count(maflat[idxFirst:idxNext]) == 0:                                 # all elements are masked
##            print "Warning: cannot do linear interpolation, all values are masked"
##        elif idx == idxFirst:                                                        # the first element is masked
##            maflat[idx] = maflat[idx+1:idxNext].compressed()[0]
##        elif idx == idxNext-1:                                                      # the last element is masked
##            maflat[idx] = maflat[idxFirst:idx].compressed()[-1]
##        else:
##            maflat[idx] = MA.average([maflat.take(range(idx-1,idxFirst-1,-1)).compressed(), maflat[idx+1:idxNext].compressed(), MA.masked])
####            maflat[idx] = (maflat[idxFirst:idx].compressed()[-1] + maflat[idx+1:idxNext].compressed()[0]) / 2.
##    return MA.reshape(maflat, ma.shape)
##
##def fillLinear(ma, axis=0):
##    """Linear interpolation of missing values;
##    """
##    pass
##
##
##si = scipy.interpolate.interp1d(Numeric.array([[0,1,3,4],[0,1,3,4]]), Numeric.array([[0.1,1.1,3.2,4.4],[0,1,3,4]]))
##si(2)

###################################################################################
## swapaxes for masked array
###################################################################################

def swapaxesMA(ma, axis1, axis2):
    if axis1 < axis2:
        m,M = axis1, axis2
    elif axis1 > axis2:
        m,M = axis2, axis1
    elif 0 <= axis1 == axis2 < len(ma.shape):
        return ma
    else:
        raise ValueError("Bad axis1 or axis2 argument to swapaxes.")
    axList = range(m) + [M] + range(m+1,M) + [m] + range(M+1,len(ma.shape))
    return MA.transpose(ma, axList)


###################################################################################
## compress MA to Numeric and return Numeric + indices of non-masked places
###################################################################################

def compressIndices(ma):
    """Returns 1D compressed Numeric array and the indices of the non-masked places.
    usage:  nu,ind = compressIndices(ma)
            nu = Numeric.elementwise_function(nu)
            ma = MA.put(ma, ind, nu)
    """
    ma = MA.asarray(ma)
    nonMaskedInd = Numeric.compress(1-Numeric.ravel(MA.getmaskarray(ma)), Numeric.arange(Numeric.multiply.reduce(ma.shape)))
    return MA.filled(ma.compressed()), nonMaskedInd


###################################################################################
## compact print of Numeric or MA
###################################################################################

def aprint(a,decimals=2):
    """prints array / masked array of floats"""
    if isinstance(a, Numeric.ArrayType):
        print Numeric.around(a.astype(Numeric.Float0),decimals)
    elif type(a) == MA.MaskedArray:
        print MA.around(a.astype(Numeric.Float0),decimals)


###################################################################################
## min / max along the given axis
## SLOW DUE TO SORT
###################################################################################

def minMA(m,axis=0):
    """slow: remove sorting"""
    m = MA.asarray(m, Numeric.Float)
    transList = [axis] + range(0,axis) + range(axis+1,MA.rank(m))
    m = MA.transpose(m, transList)    # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return MA.sort(m, 0, fill_value=1e20)[0]
    
def maxMA(m,axis=0):
    """slow: remove sorting"""
    m = MA.asarray(m, Numeric.Float)
    transList = [axis] + range(0,axis) + range(axis+1,MA.rank(m))
    m = MA.transpose(m, transList)    # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return MA.sort(m, 0, fill_value=-1e20)[-1]


###################################################################################
## median absolute deviation
###################################################################################

def mad(m,axis=0):
    """Returns Median Absolute Deviation of the given array along the given axis.
    """
    m = Numeric.asarray(m)
    mx = Numeric.asarray(median(m,axis),Numeric.Float)
    xt = Numeric.transpose(m, [axis]+range(axis)+range(axis+1,Numeric.rank(m))) # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return MLab.median(Numeric.absolute(xt-mx))

def madMA(m,axis=0):
    """Returns Median Absolute Deviation of the given masked array along the given axis.
    """
    m = MA.asarray(m)
    mx = MA.asarray(medianMA(m,axis),Numeric.Float)
    xt = MA.transpose(m, [axis]+range(axis)+range(axis+1,MA.rank(m)))   # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return medianMA(MA.absolute(xt-mx))



if __name__ == "__main__":
##    n2 = numpy.asarray([[1,2],[3,4],[5,6]], numpy.float)    #shape (3,2)
##    n3 = numpy.asarray([n2, 2*n2, 3*n2, 4*n2])              #shape (4,3,2)
##
##    m2 = numpy.ma.asarray(n2)
##    m2[0,0] = numpy.ma.masked
##    m2[1,1] = numpy.ma.masked
##
##    m3 = numpy.ma.asarray(n3)
##    m3[0,0,0] = 10000
##    m3[1,1,1] = 10000
##    m3[2,2,0] = 10000
##    m3[3,0,1] = 10000
##    m3[0,0,0] = numpy.ma.masked
##    m3[1,1,1] = numpy.ma.masked
##    m3[2,2,0] = numpy.ma.masked
##    m3[3,0,1] = numpy.ma.masked
##
##    mm3 = numpy.ma.ones((4,3,2), numpy.float) * numpy.ma.masked
##    mm3[0,0,0] = 33
##    mm3[1,1,1] = 333
##    mm3[2,2,0] = 3333
##    mm3[3,0,1] = 33333
##
##    # compare to Numeric
##    import Numeric, MA
##    n2n = Numeric.asarray([[1,2],[3,4],[5,6]], Numeric.Float)    #shape (3,2)
##    n3n = Numeric.asarray([n2n, 2*n2n, 3*n2n, 4*n2n])              #shape (4,3,2)
##
##    m2n = MA.asarray(n2n)
##    m2n[0,0] = MA.masked
##    m2n[1,1] = MA.masked
##
##    m3n = MA.asarray(n3n)
##    m3n[0,0,0] = 10000
##    m3n[1,1,1] = 10000
##    m3n[2,2,0] = 10000
##    m3n[3,0,1] = 10000
##    m3n[0,0,0] = MA.masked
##    m3n[1,1,1] = MA.masked
##    m3n[2,2,0] = MA.masked
##    m3n[3,0,1] = MA.masked
##
##    mm3n = MA.ones((4,3,2), Numeric.Float) * MA.masked
##    mm3n[0,0,0] = 33
##    mm3n[1,1,1] = 333
##    mm3n[2,2,0] = 3333
##    mm3n[3,0,1] = 33333

    # argmaxMA    
    m2 = numpy.ma.asarray([[1,2,3],[4,5,6]])*numpy.ma.masked
    m2[1,2] = 10    
##    print "orig, axis None\t", numpy.ma.argmax(m2, None)
##    print "extn, axis None\t", argmaxMA(m2, None)
##    print "orig, axis 0\t", numpy.ma.argmax(m2,0)
##    print "extn, axis 0\t", argmaxMA(m2,0)
##    print "orig, axis 1\t", numpy.ma.argmax(m2,1)
##    print "extn, axis 1\t", argmaxMA(m2,1)
    ma = numpy.ma.argmax(m2,1)
    mb = argmaxMA(m2,1)
