import types, string
import Numeric, MA
import statc
    

def kNNimputeMA(arr2d, K=20, callback=None):
    """Returns a new 2D MA.array with missing values imputed from K nearest neighbours.
    Find K rows (axis 0) with the most similar values in axis 1; similarity measure: weighted Euclidean distance.
    Imputed value = weighted average of the corresponding values of K nearest neighbours,
    where weights equal to tricubic distribution of distances to all rows.
    """
    arr2d = MA.asarray(arr2d)
    assert len(arr2d.shape) == 2, "2D array expected"
    # make a copy for imputation
    aImp = MA.array(arr2d)
    # impute only values in the columns that have more than K non-masked values
    countByColumns = MA.count(arr2d, axis=0)
    columnCond = Numeric.greater(countByColumns, K)
    columnInd = Numeric.arange(arr2d.shape[1])
    columnIndCompress = Numeric.compress(columnCond, columnInd) # column indices where the imputation in possible (there are at least K+1 known values)
    # impute the rows where 1 < #known_values < len(row), i.e. exclude the rows with 0, 1 and all values
    countByRows = MA.count(MA.take(arr2d, columnIndCompress, axis=1), axis=1)
    rowCond = Numeric.logical_and(Numeric.greater(countByRows, 1), Numeric.less(countByRows, len(columnIndCompress)))
    # start imputation rowwise
    for rowIdx in Numeric.compress(rowCond, Numeric.arange(arr2d.shape[0])):
        rowResized = MA.resize(arr2d[rowIdx], arr2d.shape)
        diff = arr2d - rowResized
        numComp = MA.count(diff, axis=1)
        distances = MA.sqrt(MA.add.reduce((diff)**2, 1) / numComp)
        indSorted = MA.argsort(distances, fill_value=1e20)[1:]  # nearest neighbours row indices (without the current row index)
        distSort = MA.take(distances, indSorted)
        distSortCompress = MA.compress(-1*(distSort.mask()-1), distSort).filled()
        weightsSortCompress = Numeric.power(1-Numeric.power(distSortCompress/distSortCompress[-1],3),3) # tricubic distribution of all weights
        avrg = MA.average(MA.take(arr2d, indSorted[0:K], axis=0), weights=weightsSortCompress[0:K]) # average from K neighbours
        columnInd4currentRow = Numeric.compress(Numeric.logical_and(arr2d[rowIdx].mask(), columnCond), columnInd)
        for colIdx in columnInd4currentRow:
            aImp[rowIdx,colIdx] = avrg[colIdx]
        if callback:
            callback()
    return aImp


def loessMA(m, windowSize, axis=0, callback=None):
    """Returns a new array with values at the given axis smoothed by loess; the masked values are approximated by loess.
    Assumes equidistant spacing of points on the given axis.
    """
    assert 0 < windowSize <= m.shape[axis], "0 < windowSize_% <=1 OR windowSize in range(1,m.shape[axis]+1)"
    m = MA.asarray(m)
    shp_other = list(m.shape)
    shp_other.pop(axis)
    # get a transposed and reshaped mask and data from m; if m.mask() == None, construct a new array of zeros
    mask = Numeric.reshape(Numeric.transpose(MA.getmaskarray(m), [axis] + range(0,axis) + range(axis+1,len(m.shape))), (m.shape[axis], Numeric.multiply.reduce(shp_other)))
    data = Numeric.reshape(Numeric.transpose(m.raw_data(), [axis] + range(0,axis) + range(axis+1,len(m.shape))), (m.shape[axis], Numeric.multiply.reduce(shp_other)))
    maskInv = -1*(mask-1)
    xall = Numeric.arange(data.shape[0])
    for ii in range(data.shape[1]):
        data[:,ii] = Numeric.array(statc.loess(zip(Numeric.compress(maskInv[:,ii], xall).tolist(), Numeric.compress(maskInv[:,ii], data[:,ii]).tolist()), xall.tolist(), windowSize))[:,1]
        if callback:
            callback()
    return Numeric.transpose(Numeric.reshape(data, [m.shape[axis]] + shp_other), [axis] + range(0,axis) + range(axis+1,len(m.shape)))


if __name__=="__main__":
    import Dicty
    DN = Dicty.DData.DData_Nancy()
    y2 = DN.getRaw2d("yakA1.2.raw")
    y2i = kNNimputeMA(y2)
