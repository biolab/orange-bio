import math
import Numeric
import MA
import MLab, LinearAlgebra

import numarray as NA
import numarray.ma as NM
import numarray.linear_algebra as LA


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
            nu = Numeric.f_elementwise(nu)
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
    """MA -> numarray.ma"""
    m = MA.asarray(m)
    return NM.array(m.raw_data(), mask=m.raw_mask())


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
    return MA.array(m.raw_data(), mask=m.raw_mask())


###################################################################################
## min / max along the given axis
## SLOW DUE TO SORT
###################################################################################

def minMA(m,axis=0):
    """slow: remove sorting"""
    m = MA.asarray(m, MA.Float)
    transList = [axis] + range(0,axis) + range(axis+1,MA.rank(m))
    m = MA.transpose(m, transList)    # do not use swapaxes
    return MA.sort(m, 0, fill_value=1e20)[0]
    
def maxMA(m,axis=0):
    """slow: remove sorting"""
    m = MA.asarray(m, MA.Float)
    transList = [axis] + range(0,axis) + range(axis+1,MA.rank(m))
    m = MA.transpose(m, transList)    # do not use swapaxes
    return MA.sort(m, 0, fill_value=-1e20)[-1]

###################################################################################
## median
###################################################################################

def median(m,axis=0):
    m = Numeric.asarray(m)
    return MLab.median(Numeric.transpose(m, [axis]+range(axis)+range(axis+1,Numeric.rank(m))))  # do not use swapaxes


def medianMA(m,axis=0):
    """Returns median of the given masked array along the given axis."""
    return _nm2ma(_medianNM(_ma2nm(m),axis))


def _medianNA(m,axis=0):
    m = NA.asarray(m)
    return LA.mlab.median(NA.transpose(m, [axis]+range(axis)+range(axis+1,m.rank)))   # do not use swapaxes

def _medianNM(m,axis=0):
    """Returns median of the given masked numarray along the given axis."""
    m = NM.asarray(m, NM.Float)
    if len(m.shape) == 0:
        return m[0]
    elif len(m.shape) == 1:
        c = NM.count(m)
        mIndSort = NM.argsort(m,fill_value=1e20)
        return (m[mIndSort[(c-1)/2]] + m[mIndSort[c/2]]) / 2.
    else:
        # count the number of nonmasked indices along the given axis
        c = NM.count(m, axis=axis)
        # prepare indices for other axis (except for the given)
        cind = NA.indices(NA.shape(c))
        # indices of m sorted by increasing value at axis 0; masked values placed at the end;
        # use _argsortSwapNM() to get the shape same as m
        mtIndSort = _argsortSwapNM(m,axis,fill_value=1e20)
        # get indices to address median elements from m
        # tuple(takeInd1): all lists in the tuple must be of the same shape,
        #   the result mtIndSort[tuple(takeInd1)] is also of that shape
        takeInd1 = cind.tolist()
        takeInd1.insert(axis,((c-1)/2).tolist())
        medInd1 = mtIndSort[tuple(takeInd1)]
        takeInd2 = cind.tolist()
        takeInd2.insert(axis,(c/2).tolist())
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
        return (med1+med2)/2.

def _argsortSwapNM(m,axis=0,fill_value=None):
    """Returns the indices along the given axis sorted by increasing value of the given masked numarray
    e.g. m=[[.5,.2],[.3,.6]], _argsortSwapNM(m,0)->[[1,0],[0,1]]"""
    m = NM.asarray(m)
    return NA.swapaxes(NA.asarray(NM.argsort(m, axis, fill_value)), -1, axis)


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
    xt = Numeric.transpose(m, [axis]+range(axis)+range(axis+1,Numeric.rank(m))) # do not use swapaxes
    return MLab.median(Numeric.absolute(xt-mx))

def madMA(m,axis=0):
    """Returns Median Absolute Deviation of the given masked array along the given axis.
    """
    m = MA.asarray(m)
    mx = MA.asarray(medianMA(m,axis),MA.Float)
    xt = MA.transpose(m, [axis]+range(axis)+range(axis+1,MA.rank(m)))   # do not use swapaxes
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


###################################################################################
## binomial factors
###################################################################################

def binomial(n,r):
    """Returns array of elementwise binomial factors: n[i]! / r[i]!(n[i]-r[i])!
    type(n) & type(r) == int | Numeric.ArrayType of int
    compute in floats to avoid overflow; max. n = 1029
    """
##    n = Numeric.asarray(n).astype(Numeric.Int)
##    r = Numeric.asarray(r).astype(Numeric.Int)
    n = Numeric.asarray(n)
    r = Numeric.asarray(r)
    assert Numeric.shape(n) == Numeric.shape(r), "binomial(n,r): n and r of different shape"
    bin = Numeric.zeros(Numeric.shape(n), Numeric.Float)
    for idx in xrange(len(n.flat)):
        ni1 = n.flat[idx] + 1
        b = Numeric.ones((ni1), Numeric.Float)
        for i in xrange(1,ni1):
            for j in xrange(i-1,0,-1):
                b[j] += b[j-1]
        bin.flat[idx] = b[r.flat[idx]]
    return bin

def __binomial_single(n,r):
    """
    DEPRICATED, see binomial(n,r)
    returns n! / r!(n-r)!
    """
    b = Numeric.ones((n+1), Numeric.Float)
    for i in xrange(1,n+1):
        for j in xrange(i-1,0,-1):
            b[j] += b[j-1]
    return b[r]






if __name__ == "__main__":

    m1 = NM.array([4,3,6,8,54,1], mask=[1,1,0,1,1,1])
    m2a = NM.array([[1,2,3],[4,5,6],[7,8,9]], NM.Float, mask=[[0,1,0],[0,1,0],[0,0,1]])
    n2 = NA.array([[3,1,2],[4,3,6],[1,4,7],[4,3,6]])
    n3 = NA.array([n2,n2])
    m2 = NM.array([[5,7,3,3],[6,4,5,1],[2,6,1,4],[3,4,5,1],[3,4,5,1],[6,4,2,3]], NM.Int, mask=[[0,1,0,0],[0,1,0,0],[0,0,0,1],[0,0,0,0],[0,0,0,1],[1,1,1,1]])
    m2s = NM.array([[5,7,3,3],[6,4,5,1],[2,6,1,4],[3,4,5,1],[3,4,5,1],[6,4,2,3],[6,4,2,3],[6,4,2,3],[6,4,2,3],[6,4,2,3]], NM.Int, mask=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])
    m2u = NM.array([[5,7,3,3],[6,4,5,1],[2,6,1,4],[3,4,5,1],[3,4,5,1],[6,4,2,3]], NM.Int, mask=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    m3 = NM.array([m2.raw_data(),10*m2.raw_data()],NM.Int, mask=[m2.raw_mask(),m2.raw_mask()], fill_value=-1)

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
    