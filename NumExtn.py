import math
import Numeric
import MA
import MLab, LinearAlgebra

import numarray as NA
import numarray.ma as NM
import numarray.linear_algebra as LA
import numpy


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

from Permutation import *

def permutInverse(n):
    """Returns inverse permutation given integers in range(len(n)),
    such that permitInverse(permutInverse(range(4)))==range(4).
    """
    n = Numeric.asarray(n)
    pInv = Numeric.argsort(n)
    assert Numeric.equal(n, Numeric.argsort(pInv)), "Inverse not successful; input should be permutation of range(len(input))."
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
    r = MA.zeros(m.shape[0], MA.Float)
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

    m2d = MA.ravel(MA.zeros((m2dShape0, m2dShape0), m1d.typecode()))
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
    return MA.take(MA.ravel(m2d), takeInd)


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
##            maflat[idx] = MA.average([MA.take(maflat,range(idx-1,idxFirst-1,-1)).compressed(), maflat[idx+1:idxNext].compressed(), MA.masked])
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
    if type(a) == Numeric.ArrayType:
        print Numeric.around(a.astype(Numeric.Float0),decimals)
    elif type(a) == MA.MaskedArray:
        print MA.around(a.astype(MA.Float0),decimals)


###################################################################################
## Numeric -> numarray
## MA -> numarray.ma
###################################################################################

def _nu2na(a):
    """Numeric -> numarray"""
    a = Numeric.asarray(a)
    return NA.array(a)

def _ma2nm(m):
    """MA -> numarray.ma
    BUG: crashes with the following array:
                      array(data =
                 [[1,2,3,]
                 [0,0,0,]
                 [7,8,9,]],
                      mask =
                 [[0,0,0,]
                 [1,1,1,]
                 [0,0,0,]],
                      fill_value=[0,])
"""
    m = MA.asarray(m)
    m = NM.array(m.raw_data(), mask=m.raw_mask())
    m.tolist()  # execute this line in order to fix the discribed bug
    return m


###################################################################################
## numarray -> Numeric
## numarray.ma -> MA
###################################################################################

def _na2nu(a):
    """numarray -> Numeric"""
    a = NA.asarray(a)
    return Numeric.array(a)

def _nm2ma(m):
    """numarray.ma -> MA"""
    m = NM.asarray(m)
    return MA.array(m.raw_data(), m.typecode(), mask=m.raw_mask())


###################################################################################
## min / max along the given axis
## SLOW DUE TO SORT
###################################################################################

def minMA(m,axis=0):
    """slow: remove sorting"""
    m = MA.asarray(m, MA.Float)
    transList = [axis] + range(0,axis) + range(axis+1,MA.rank(m))
    m = MA.transpose(m, transList)    # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return MA.sort(m, 0, fill_value=1e20)[0]
    
def maxMA(m,axis=0):
    """slow: remove sorting"""
    m = MA.asarray(m, MA.Float)
    transList = [axis] + range(0,axis) + range(axis+1,MA.rank(m))
    m = MA.transpose(m, transList)    # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return MA.sort(m, 0, fill_value=-1e20)[-1]

###################################################################################
## median
###################################################################################

def median(m,axis=0):
    m = Numeric.asarray(m)
    return MLab.median(Numeric.transpose(m, [axis]+range(axis)+range(axis+1,Numeric.rank(m))))  # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)


def medianMA(m,axis=0):
    """Returns median of the given masked array along the given axis."""
    return _nm2ma(_percentilesNM(_ma2nm(m),0.5,axis))


def _medianNA(m,axis=0):
    m = NA.asarray(m)
    return LA.mlab.median(NA.transpose(m, [axis]+range(axis)+range(axis+1,m.rank)))   # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)

##def _medianNM(m,axis=0):
##    """Returns median of the given masked numarray along the given axis."""
##    m = NM.asarray(m, NM.Float)
##    if len(m.shape) == 0:
##        return m[0]
##    elif len(m.shape) == 1:
##        c = NM.count(m)
##        mIndSort = NM.argsort(m,fill_value=1e20)
##        return (m[mIndSort[(c-1)/2]] + m[mIndSort[c/2]]) / 2.
##    else:
##        # count the number of nonmasked indices along the given axis
##        c = NM.count(m, axis=axis)
##        # prepare indices for other axis (except for the given)
##        cind = NA.indices(NA.shape(c))
##        # indices of m sorted by increasing value at axis 0; masked values placed at the end;
##        # use _argsortSwapNM() to get the shape same as m
##        mtIndSort = _argsortSwapNM(m,axis,fill_value=1e20)
##        # get indices to address median elements from m
##        # tuple(takeInd1): all lists in the tuple must be of the same shape,
##        #   the result mtIndSort[tuple(takeInd1)] is also of that shape
##        takeInd1 = cind.tolist()
##        takeInd1.insert(axis,((c-1)/2).tolist())
##        medInd1 = mtIndSort[tuple(takeInd1)]
##        takeInd2 = cind.tolist()
##        takeInd2.insert(axis,(c/2).tolist())
##        medInd2 = mtIndSort[tuple(takeInd2)]
##        # get both median elements; if c[i] odd: med1[i] == med2[i]
##        takeMed1 = cind.tolist()
##        takeMed1.insert(axis,medInd1.tolist())
##        med1 = m[tuple(takeMed1)]
##        takeMed2 = cind.tolist()
##        takeMed2.insert(axis,medInd2.tolist())
##        med2 = m[tuple(takeMed2)]
####        if __name__=="__main__":
####            print "m\n",m.filled(-1)
####            print "c", c
####            print "[(c-1)/2,c/2]", [(c-1)/2,c/2]
####            print "cind\n",cind
####            print "mtIndSort\n",mtIndSort
####            print "medInd1\n",medInd1
####            print "medInd2\n",medInd2
####            print "med1\n",med1
####            print "med2\n",med2
##        return (med1+med2)/2.


def percentilesMA(m,perc,axis=0):
    """Returns the percentiles of the given masked array along the given axis."""
    return _nm2ma(_percentilesNM(_ma2nm(m),perc,axis))


def _percentilesNM(m,perc,axis=0):
    """Returns the percentiles of the given masked numarray along the given axis.
    """
    assert 0 < perc < 1
    m = NM.asarray(m, NM.Float)
    if len(m.shape) == 0 or (len(m.shape)==1 and m.shape[0]==1):
        return m[0]
    elif len(m.shape) == 1:
        # returns (1,) NM array or NM.masked
        mCount = NM.count(m)
        if mCount == 0:
            return NM.masked
        elif mCount == 1:
            return m.compressed()[0]
        else:
            _k = float(perc)*(mCount-1)
            k = int(_k)
            d = _k - k
            mIndSort = NM.argsort(m,fill_value=1e20)
            return m[mIndSort[k]] + d*(m[mIndSort[k+1]]-m[mIndSort[k]])
    else:
        # count the number of nonmasked indices along the given axis
        _k = float(perc) * (NM.count(m, axis=axis)-1)
        k = _k.astype(NM.Int)
        d = _k - k
        # prepare indices for other axis (except for the given)
        cind = NA.indices(NA.shape(k))
        # indices of m sorted by increasing value at axis 0; masked values placed at the end;
##        2006-08-25: PJ: _argsortSwapNM replaced by NM.argsort
##        # use _argsortSwapNM() to get the shape same as m
##        mtIndSort = _argsortSwapNM(m,axis,fill_value=1e20)
        mtIndSort = NM.argsort(m,axis,fill_value=1e20)
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
##            print "m\n",m.filled(-1)
##            print "c", c
##            print "[(c-1)/2,c/2]", [(c-1)/2,c/2]
##            print "cind\n",cind
##            print "mtIndSort\n",mtIndSort
##            print "medInd1\n",medInd1
##            print "medInd2\n",medInd2
##            print "med1\n",med1
##            print "med2\n",med2
        return med1 + d*(med2-med1).filled(0)


##def _argsortSwapNM(m,axis=0,fill_value=None):
##    """Returns the indices along the given axis sorted by increasing value of the given masked numarray
##    e.g. m=[[.5,.2],[.3,.6]], _argsortSwapNM(m,0)->[[1,0],[0,1]]
##
##    Used only for numarray ver. 0.9 to correct a BUG in function argsort(...)
##    """
##    m = NM.asarray(m)
##    return NA.swapaxes(NA.asarray(NM.argsort(m, axis, fill_value)), -1, axis)


###################################################################################
## std. deviation MA
###################################################################################

def stdMA(m,axis=0):
    m = MA.asarray(m)
    m = MA.transpose(m, [axis] + range(0,axis) + range(axis+1,MA.rank(m)))  # do not use swapaxes, use (axis, 0...axis-1, axis+1...rank)
    return MA.sqrt(MA.add.reduce((MA.subtract(m,MA.average(m,0)))**2,0)/(m.shape[0]-1))


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
    mx = MA.asarray(medianMA(m,axis),MA.Float)
    xt = MA.transpose(m, [axis]+range(axis)+range(axis+1,MA.rank(m)))   # do not use swapaxes: (0,1,2) -swap-> (2,1,0); (0,1,2) -transpose-> (2,0,1)
    return medianMA(MA.absolute(xt-mx))

def _madNA(m,axis=0):
    """Returns Median Absolute Deviation of the given numarray along the given axis.
    """
    m = NA.asarray(m)
    mx = NA.asarray(_medianNA(m,axis),NA.Float)
    xt = NA.transpose(m, [axis]+range(axis)+range(axis+1,m.rank))
    return LA.mlab.median(NA.absolute(xt-mx))

def _madNM(m,axis=0):
    """Returns Median Absolute Deviation of the given masked numarray along the given axis.
    """
    m = NM.asarray(m, NM.Float)
    mx = _medianNM(m, axis)
    # BUG WARNING: transpose returns ObjectArray, use NM.array(..., NM.Float)
    xt = NM.array(NM.transpose(m, [axis]+range(axis)+range(axis+1,NM.rank(m))), NM.Float)
    return _medianNM(NM.absolute(xt-mx),0)









if __name__ == "__main__":
    pass

##    m1 = NM.array([4,3,6,8,54,1], mask=[1,1,0,1,1,1])
##    m2a = NM.array([[1,2,3],[4,5,6],[7,8,9]], NM.Float, mask=[[0,1,0],[0,1,0],[0,0,1]])
##    n2 = NA.array([[3,1,2],[4,3,6],[1,4,7],[4,3,6]])
##    n3 = NA.array([n2,n2])
##    m2 = NM.array([[5,7,3,3],[6,4,5,1],[2,6,1,4],[3,4,5,1],[3,4,5,1],[6,4,2,3]], NM.Int, mask=[[0,1,0,0],[0,1,0,0],[0,0,0,1],[0,0,0,0],[0,0,0,1],[1,1,1,1]])
##    m2s = NM.array([[5,7,3,3],[6,4,5,1],[2,6,1,4],[3,4,5,1],[3,4,5,1],[6,4,2,3],[6,4,2,3],[6,4,2,3],[6,4,2,3],[6,4,2,3]], NM.Int, mask=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])
##    m2u = NM.array([[5,7,3,3],[6,4,5,1],[2,6,1,4],[3,4,5,1],[3,4,5,1],[6,4,2,3]], NM.Int, mask=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
##    m3 = NM.array([m2.raw_data(),10*m2.raw_data()],NM.Int, mask=[m2.raw_mask(),m2.raw_mask()], fill_value=-1)

##    def test3D(m3):
##        print "median\n",_medianNM(m3,0).filled(-1)
##        print "......."
##        for ii in range(m3.shape[0]):
##            print m3[ii,:,:].filled(-1)
##            print
##        print "median\n",_medianNM(m3,1).filled(-1)
##        print "......."
##        for ii in range(m3.shape[1]):
##            print m3[:,ii,:].filled(-1)
##            print
##        print "median\n",_medianNM(m3,2).filled(-1)
##        print "......."
##        for ii in range(m3.shape[2]):
##            print m3[:,:,ii].filled(-1)
##            print

    # median    
##    m1med = _medianNM(m1)
##    m2med = _medianNM(m2,1)
##    m3med = _medianNM(m3,2)

    # MAD: must be equal!
##    mad(m2u)
##    mad(m2s)

    # compare Numeric, MA, NA, NM
##    import RandomArray
##    lst1 = (RandomArray.random((7)) * 10).tolist()
##    lst2 = (RandomArray.random((4,7)) * 10).tolist()
##    lst3 = (RandomArray.random((3,4,7)) * 10).tolist()
##
##    def getAllArrays(lst1):
##        # Numeric, MA
##        nu1 = Numeric.array(lst1)
##        ma1 = MA.array(lst1)
##        ma1m = MA.masked_less(lst1, 5)
##        # numarray, ma
##        na1 = NA.array(lst1)
##        nm1 = NM.array(lst1)
##        nm1m = NM.masked_less(lst1, 5)
##        return nu1, ma1, ma1m, na1, nm1, nm1m
##
##    nu1, ma1, ma1m, na1, nm1, nm1m = getAllArrays(lst1)
##    nu2, ma2, ma2m, na2, nm2, nm2m = getAllArrays(lst2)
##    nu3, ma3, ma3m, na3, nm3, nm3m = getAllArrays(lst3)

    # linear interpolation of missing values (fillLinear)
##    m1 = MA.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]], mask = [[1,0,0,0],[0,1,0,0],[0,0,0,1]])
##    m2 = MA.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]], mask = [[1,1,0,0],[0,1,1,0],[0,0,1,1]])
##    m1l = fillLinear(m1,1)
##    m2l = fillLinear(m2,1)

    # triangular matrix <-> 1D masked array
##    m1 = MA.arange(0, 6)
##    m2u = triangularPut(m1,1,0)
##    m2l = triangularPut(m1,0,1)
##    m2ul = triangularPut(m1,1,1)
##    m2 = triangularPut(m1,0,0)
##
##    m2uld = diagonalPut(Numeric.arange(4), m2ul)
##
##    m1u = triangularGet(m2ul)
##    m1l = triangularGet(m2ul, 0)

    # logical unary OR
##    m1 = MA.array([0,0,0,9,9,9,1,1,1], mask=[0,0,0,1,1,1,0,0,0])
##    m2 = MA.array([0,9,1,0,9,1,0,9,1], mask=[0,1,0,0,1,0,0,1,0])
##    m1uorm2 = logical_unary_or(m1,m2)
##    print m1uorm2 == MA.array([0,0,1,0,9,1,1,1,1], mask=[0,0,0,0,1,0,0,0,0])
##    m1uandm2 = logical_unary_and(m1,m2)
##    print m1uandm2 == MA.array([0,0,0,0,9,1,0,1,1], mask=[0,0,0,0,1,0,0,0,0])

    # apply binary
##    m1 = MA.array([0,0,0,9,9,9,1,1,1], mask=[0,0,0,1,1,1,0,0,0])
##    m2 = MA.array([0,9,1,0,9,1,0,9,1], mask=[0,1,0,0,1,0,0,1,0])
##    m3 = MA.array([1,1,0,0,1,1,0,0,1], mask=[0,1,0,1,0,1,0,1,0])
##    ors = apply_binary(logical_unary_or, m1,m2,m3)
##    ands = apply_binary(logical_unary_and, m1,m2,m3)

    
    # DEBUG: _ma2nm
##    d3 = MA.array([[1,2,3],[4,5,6],[7,8,9]])
##    d2 = d3 / [1.,0.,1.]
##    d1 = d3 / [0.,0.,1.]
##    d0 = d3 / [0.,0.,0.]
##    print "d3", _ma2nm(d3).tolist()
##    print "d2", _ma2nm(d2).tolist()
##    print "d1", _ma2nm(d1).tolist()
##    print "d0", _ma2nm(d0).tolist()
##    print _percentilesNM(_ma2nm(d3), 0.5, 0).tolist()
##    print _percentilesNM(_ma2nm(d2), 0.5, 0).tolist()
##    print _percentilesNM(_ma2nm(d1), 0.5, 0).tolist()
##    print _percentilesNM(_ma2nm(d0), 0.5, 0).tolist()
    