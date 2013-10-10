## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA
import statc
    

def kNNimputeMA(arr2d, K=20, callback=None):
    """Returns a new 2D MA.array with missing values imputed from K nearest neighbours.
    Find K rows (axis 0) with the most similar values where similarity measure corresponds to weighted Euclidean distance.
    Imputed value = weighted average of the corresponding values of K nearest neighbours,
    where weights equal to tricubic distribution of distances to all rows.
    Impute missing rows by average over all rows.
    Version: 30.8.2005
    """
    arr2d = MA.asarray(arr2d)
    assert len(arr2d.shape) == 2, "2D array expected"
    # make a copy for imputation
    aImp2 = MA.array(arr2d)
    # leave out columns with 0 known values (columnInd: non-zero columns)
    columnCond = Numeric.greater(MA.count(arr2d, axis=0), 0)
    columnIndAll = Numeric.arange(arr2d.shape[1])
    columnInd = Numeric.compress(columnCond, columnIndAll)
    # impute the rows where 0 < #known_values < #non_zero_columns, i.e. exclude the rows with 0 and all (non-zero-column) values
    countByRows = MA.count(arr2d, axis=1)
    for rowIdx in Numeric.compress(Numeric.logical_and(Numeric.greater(countByRows, 0), Numeric.less(countByRows, columnInd.shape[0])), Numeric.arange(arr2d.shape[0])):
        rowResized = MA.resize(arr2d[rowIdx], arr2d.shape)
        diff = arr2d - rowResized
        distances = MA.sqrt(MA.add.reduce((diff)**2, 1) / MA.count(diff, axis=1))
        # nearest neighbours row indices (without the current row index)
        indSorted = MA.argsort(distances)[1:]
        distSorted = distances.take(indSorted)
        # number of distances different from MA.masked
        numNonMasked = distSorted.shape[0] - Numeric.add.reduce(Numeric.asarray(MA.getmaskarray(distSorted), Numeric.Int))
        # number of distances to account for (K or less)
        if numNonMasked > 1:
            weightsSorted = MA.power(1-MA.power(distSorted/distSorted[numNonMasked-1],3),3) # tricubic distribution of all weights
        else:
            weightsSorted = Numeric.ones(distSorted.shape[0])
        # compute average for each column separately in order to account for K non-masked values
        colInd4CurrRow = Numeric.compress(Numeric.logical_and(MA.getmaskarray(arr2d[rowIdx]), columnCond), columnIndAll)
        for colIdx in colInd4CurrRow:
            # column values sorted by distances
            columnVals = arr2d[:,colIdx].take(indSorted)
            # take only those weights where columnVals does not equal MA.masked
            weightsSortedCompressed = MA.compress(1-MA.getmaskarray(columnVals), weightsSorted)
            # impute from K (or possibly less) values
            aImp2[rowIdx,colIdx] = MA.average(columnVals.compressed()[:K], weights=weightsSortedCompressed[:K])
        if callback:
            callback()
    # impute the unknown rows with average profile
    avrgRow = MA.average(arr2d, 0)
    for rowIdx in Numeric.compress(Numeric.equal(countByRows, 0), Numeric.arange(arr2d.shape[0])):
        aImp2[rowIdx] = avrgRow
        if callback:
            callback()
    return aImp2


def loessMA(m, windowSize, axis=0, approxMasked=True, verbose=False, callback=None):
    """Returns a new array with values at the given axis smoothed by loess;
    if approxMasked==True: the masked values are approximated by loess;
    assumes equidistant spacing of points on the given axis.
    """
    assert 0 < windowSize <= m.shape[axis]+0.1, "0 < windowSize[%s] <= 1 OR windowSize in range(1.1,m.shape[axis]+1) expected, got %f" % ("%", windowSize)
    m = MA.asarray(m)
    if m.dtype.char <> Numeric.Float:
        m = m.astype(Numeric.Float)
    shp_other = list(m.shape)
    shp_other.pop(axis)
    # get a transposed and reshaped mask and data from m; if m.mask() == None, construct a new array of zeros
    mask = Numeric.reshape(Numeric.transpose(MA.getmaskarray(m), [axis] + range(0,axis) + range(axis+1,len(m.shape))), (m.shape[axis], Numeric.multiply.reduce(shp_other)))
    data = MA.reshape(MA.transpose(m, [axis] + range(0,axis) + range(axis+1,len(m.shape))), (m.shape[axis], Numeric.multiply.reduce(shp_other)))
    maskInv = -1*(mask-1)
    xall = Numeric.arange(data.shape[0])
    xallList = xall.tolist()
    for ii in Numeric.compress(Numeric.add.reduce(maskInv,0) > 1, range(data.shape[1])):    # run loess if the profile contains more than 2 values
        try:
            data[:,ii] = MA.array(statc.loess(zip(MA.compress(maskInv[:,ii], xall).tolist(), MA.compress(maskInv[:,ii], data[:,ii]).tolist()), xallList, windowSize))[:,1]
        except:
            if verbose:
                print "Warning: loessMA: could not loess axis %i index %i" % (axis, ii)
        if callback:
            callback()
    if not approxMasked:
        data = MA.array(data, mask=mask)
    return MA.transpose(MA.reshape(data, [m.shape[axis]] + shp_other), [axis] + range(0,axis) + range(axis+1,len(m.shape)))



if __name__=="__main__":
    import numpy.oldnumeric.random_array as RandomArray, time

#   10x4 array    
##    m2 = MA.asarray(RandomArray.random((10, 4)))
##    m2[0,0] = MA.masked
##    m2[1,1] = MA.masked
##    m2[2,3] = MA.masked
##    m2[3,:2] = MA.masked
##    m2[4,2:] = MA.masked
##    m2[5,:] = MA.masked
##    m2[5,1] = 0.42423
##    m2[6,:] = MA.masked    
##    m2i3 = kNNimputeMA(m2, 3)
##    m2i30 = kNNimputeMA(m2, 30)

    # 0 column
##    m3 = MA.concatenate([m2, MA.zeros((10,1))*MA.masked], 1)
##    m3i3 = kNNimputeMA(m3, 3)

    # 3 rows only
##    m34 = MA.asarray(RandomArray.random((3,4)))
##    m34[0,:2] = MA.masked
##    m34[1,1:3] = MA.masked
##    m34[-1] = MA.masked
##    m34i = kNNimputeMA(m34,3)

    # large array
##    shp = (1000,10)
##    l1 = RandomArray.random(shp)
##    m1 = RandomArray.randint(0,2,shp)
##    r1 = MA.array(l1, mask=m1)
##    print time.ctime()
##    r1i = kNNimputeMA(r1)
##    print time.ctime()

