## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

from __future__ import absolute_import

import math
import os, os.path

import numpy
import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA
import numpy.oldnumeric.linear_algebra as LinearAlgebra

import orange

from . import numpyExtn

#########################################################################
## standardization of profiles / arrays
#########################################################################

def standardize_arrays(et, robust=True):
    """Returns a new example table where values are standardized to mean zero and unit variance columnwise.
    Option: parametric (mean/stddev) or non-parametric (median/mad).
    """
    ma = orng2ma(et)
    if robust:
        return ma2orng_keepClassMetas(scaleMad(centerMed(ma)), et)
    else:
        return ma2orng_keepClassMetas(scaleStd(centerMean(ma)), et)
    

def standardize_genes(et, robust=True):
    """Returns a new example table where values are standardized to mean zero and unit variance rowwise.
    Option: parametric (mean/stddev) or non-parametric (median/mad).
    """
    ma = orng2ma(et)
    if robust:
        return ma2orng_keepClassMetas(scaleMad(centerMed(ma,1),1), et)
    else:
        return ma2orng_keepClassMetas(scaleStd(centerMean(ma,1),1), et)


#########################################################################
## merge replicas by mean, median, min, max
#########################################################################

def merge_replicas(aETStruct, type):
    """Returns a list of tuples (strain, [avrg_orngET]) where aETStruct corresponds to a list of tuples (strain, [orngET1, orngET2, ...]);
    type = ["mean" | "median" | "min" | "max"]
    """
    shape = [0,0,0]
    et0 = aETStruct[0][1][0]                                            # the first example table
    shape[0] = len(et0)                                                 # number of examples (genes)
    shape[1] = len(et0.domain.attributes)                               # number of attributes (time points)
    mergedETStruct = []
    if type == "mean":
        merge_func = MA.average
    elif type == "median":
        merge_func = numpyExtn.medianMA
    elif type == "min":
        merge_func = numpyExtn.minMA
    elif type == "max":
        merge_func = numpyExtn.maxMA
    else:
        raise AttributeError, "type = ['mean' | 'median' | 'min' | 'max']"
    for st, etList in aETStruct:
        shape[2] = len(etList)
        ma3d = MA.zeros(shape, Numeric.Float)
        for idx, et in enumerate(etList):
            ma3d[:,:,idx] = orng2ma(et)
        mergedETStruct.append((st, [ma2orng_keepClassMetas(merge_func(ma3d, 2), etList[0])]))
    return mergedETStruct
        

###################################################################################
## HELPER FUNCTIONS
###################################################################################

def scaleStd(nm, axis=0):
    """Returns new masked numarray with values scaled by dividing by Median Absolute Difference: MAD = median(|val(i)-median(val(i))|)."""
    nm = numpyExtn.swapaxesMA(MA.asarray(nm, Numeric.Float), 0, axis)
    return numpyExtn.swapaxesMA(nm / numpyExtn.stdMA(nm,0), 0, axis)

def scaleMad(nm, axis=0):
    """Returns new masked numarray with values scaled by dividing by Median Absolute Difference: MAD = median(|val(i)-median(val(i))|)."""
    nm = numpyExtn.swapaxesMA(MA.asarray(nm, Numeric.Float), 0, axis)
    return numpyExtn.swapaxesMA(nm / numpyExtn.madMA(nm,0), 0, axis)


def centerMean(nm, axis=0):
    """Returns new masked numarray with values centered by subtracting Median."""
    nm = numpyExtn.swapaxesMA(MA.asarray(nm, Numeric.Float), 0, axis)
    return numpyExtn.swapaxesMA(nm - MA.average(nm,0), 0, axis)


def centerMed(nm, axis=0):
    """Returns new masked numarray with values centered by subtracting Median."""
    nm = numpyExtn.swapaxesMA(MA.asarray(nm, Numeric.Float), 0, axis)
    return numpyExtn.swapaxesMA(nm - numpyExtn.medianMA(nm,0), 0, axis)


def getETStruct(path):
    """Returns a list of tuples (strain, [orngET1, orngET2, ...])
    """
    etStruct = []
    strains = os.listdir(path)
    for st in strains:
        pathSt = path + "\\" + st
        replicas = os.listdir(pathSt)
        stEts = []
        for rep in replicas:
            stEts.append(orange.ExampleTable(pathSt + "\\" + rep))
        etStruct.append((st, stEts))
    return etStruct


def etStruct2ma3d(aETStruct):
    """Converts a list of tuples (strain, [orngET1, orngET2, ...]) to a 3D masked array and returns it.
    """
    shape = [0,0,0]
    et0 = aETStruct[0][1][0]                                            # the first example table
    shape[0] = len(et0)                                                 # number of examples (genes)
    shape[1] = len(et0.domain.attributes)                               # number of attributes (time points)
    shape[2] = Numeric.add.reduce(map(lambda x: len(x[1]), aETStruct))  # number of ETs (replicas over all strains)
    ma3d = MA.zeros(shape, Numeric.Float)
    k = 0
    for st, etList in aETStruct:
        for et in etList:
            ma3d[:,:,k] = orng2ma(et)
            k += 1
    return ma3d


def orng2ma(aExampleTable):
    """Converts orange.ExampleTable to MA.array based on the attribute values.
    rows correspond to examples, columns correspond to attributes, class values are left out
    missing values and attributes of types other than orange.FloatVariable are masked
    """
    vals = aExampleTable.native(0, substituteDK="?", substituteDC="?", substituteOther="?")
    ma = MA.array(vals, Numeric.PyObject)
    if aExampleTable.domain.classVar != None:
        ma = ma[:,:-1]
    mask = MA.where(MA.equal(ma, "?"), 1, 0)
    for varIdx, var in enumerate(aExampleTable.domain.attributes):
        if type(var) != orange.FloatVariable:
            mask[:,varIdx] = Numeric.ones(len(aExampleTable))
    return MA.array(MA.array(ma, Numeric.PyObject, mask=mask).filled(1e20), Numeric.Float, mask=mask)


##def ma2orng(arr2d, aDomain):
##    """Converts MA.array to orange.ExampleTable where examples correspond to rows and attributes to columns.
##    masked values converted to "?" (don't knows)
##    arr2d.shape[1] must be equal to the number of the attributes of the given domain
##    domain attributes mut be of type orange.FloatVariable
##    """
##    arr2d = MA.asarray(arr2d, Numeric.PyObject)
##    assert MA.rank(arr2d) == 2, "2d array expected"
##    assert len(aDomain.attributes) == arr2d.shape[1], "the shape of the array incompatible with the given domain"
##    if aDomain.classVar != None:
##        et = orange.ExampleTable(orange.Domain(aDomain.attributes, None), arr2d.tolist("?"))
##        return orange.ExampleTable(aDomain, et)
##    else:
##        return orange.ExampleTable(aDomain, arr2d.tolist("?"))

def ma2orng(arr2d, aDomain):
    """Converts MA.array to orange.ExampleTable where examples correspond to rows and attributes to columns.
    masked values converted to "?" (don't knows)
    arr2d.shape[1] must be equal to the number of the attributes of the given domain
    domain attributes mut be of type orange.FloatVariable
    2009-06-23: adapted to numpy
    """
    arr2d = numpy.ma.asarray(arr2d, numpy.object)
    assert numpy.ma.rank(arr2d) == 2, "2d array expected"
    assert len(aDomain.attributes) == arr2d.shape[1], "the shape of the array incompatible with the given domain"
    if aDomain.classVar != None:
        et = orange.ExampleTable(orange.Domain(aDomain.attributes, None), arr2d.tolist("?"))
        return orange.ExampleTable(aDomain, et)
    else:
        return orange.ExampleTable(aDomain, arr2d.tolist("?"))


def ma2orng_keepClassMetas(arr2d, aExampleTable):
    """Creates new example table where attribute values correspond to the given 2D array, class and meta attributes remain unchanged.
    """
    arr2d = MA.asarray(arr2d, Numeric.PyObject)
    assert MA.rank(arr2d) == 2, "2D array expected"
    assert arr2d.shape[0] == len(aExampleTable), "arr2d.shape[0] != len(aExampleTable)"
    assert arr2d.shape[1] == len(aExampleTable.domain.attributes), "arr2d.shape[1] != len(aExampleTable.domain.attributes)"
    domAtt = orange.Domain(aExampleTable.domain.attributes, None)
    if aExampleTable.domain.classVar != None:
        domClassMeta = orange.Domain([aExampleTable.domain.classVar])
    else:
        domClassMeta = orange.Domain([])
    domClassMeta.addmetas(aExampleTable.domain.getmetas())
    etAtt = orange.ExampleTable(domAtt, arr2d.tolist("?"))
    etClassMeta = orange.ExampleTable(domClassMeta, aExampleTable)
    return orange.ExampleTable(aExampleTable.domain, orange.ExampleTable([etAtt, etClassMeta]))



if __name__ == "__main__":
    
    import Meda.Preproc
    reload(Meda.Preproc)
    import Dicty
    reload(Dicty)
    import Dicty.DAnnotation

##    myPath = r"C:\Python23\Lib\site-packages\orange\OrangeWidgets\Genomics\chipdata"
##    myPath = r"C:\Documents and Settings\peterjuv\My Documents\Transcriptional Phenotype\DictyData\orngDataStruct_Nancy"
    myPath = r"C:\Documents and Settings\peterjuv\My Documents\Transcriptional Phenotype\DictyData\orngDataStruct_Nancy_yakApufA"

    def generateETStruct(path, medaData, numGenes=None):
        ddbList = Dicty.DAnnotation.getDDBList()
        if not os.path.exists(path):
            os.mkdir(path)
        medaData = Dicty.DData.DData_Nancy()    
        for st in medaData.strains:
            pathSt = path + "\\" + st
            if not os.path.exists(pathSt):
                os.mkdir(pathSt)
            for rep in medaData.strain2replicaList(st):
                ma2d = medaData.getRaw2d(rep)
                et = Meda.Preproc.ma2orng(ma2d,Meda.Preproc.getTcDomain(ma2d.shape[1], False, [], None))
                et.domain.addmeta(orange.newmetaid(), orange.StringVariable("DDB"))
                for eIdx,e in enumerate(et):
                    e["DDB"] = ddbList[eIdx]
                if numGenes:
                    orange.saveTabDelimited(pathSt + "\\" + rep + ".tab", orange.ExampleTable(et[:numGenes]))
                else:
                    orange.saveTabDelimited(pathSt + "\\" + rep + ".tab", et)

    # generate orange example table structure
##    generateETStruct(myPath, Dicty.DData.BR_ACS())

    # load data
##    etStruct = getETStruct(myPath)
##    DN = Dicty.DData.DData_Nancy()
##    rglc = DN.getReplicaGroupListCopy()
##    DNpufAyakA = DN.getSubData([rglc[DN.getIdxStrain("pufA")], rglc[DN.getIdxStrain("yakA")], rglc[DN.getIdxStrain("yakApufA")]])
##    ma3d = DNpufAyakA.getRaw3dAll()

    # test orngET <-> MA
##    m2d1 = DNpufAyakA.getRaw2d(DNpufAyakA.replicas[1])
##    et = ma2orng(m2d1, Meda.Preproc.getTcDomain(13, False, [], None))
##    m2d2 = orng2ma(et)
##    print "conversion orngET <-> MA OK?", MA.add.reduce(MA.add.reduce(m2d1 - m2d2 > 0.001)) == 0
##    print "masked places OK?", MA.count(m2d1) == MA.count(m2d2)

    # test conversion to etStruct
##    print "conversion to etStruct OK?", MA.add.reduce(MA.add.reduce(MA.add.reduce((etStruct2ma3d(etStruct) - DNpufAyakA.getRaw3dAll() > 0.001)))) == 0

    # test standardization
##    etA = standardize_arrays(etStruct[0][1][1])
##    print "std.arrays OK?", MA.add.reduce(MA.add.reduce(orng2ma(etA) - Meda.Preproc.scaleMad(Meda.Preproc.centerMed(ma3d[:,:,1])) > 0.1)) == 0
##    etB = standardize_genes(etStruct[0][1][1])
##    print "std.genes OK?", MA.add.reduce(MA.add.reduce(orng2ma(etB) - Meda.Preproc.scaleMad(Meda.Preproc.centerMed(ma3d[:,:,1],1),1) > 0.1)) == 0
##    
    # test merge replicas
##    mergedETStruct_mean = merge_replicas(etStruct, "mean")
##    print "merge replicas: mean OK?", MA.add.reduce(MA.add.reduce(orng2ma(mergedETStruct_mean[0][1][0]) - DNpufAyakA.getMean2d("pufA") > 0.001)) == 0
##    mergedETStruct_median = merge_replicas(etStruct, "median")
##    mergedETStruct_min = merge_replicas(etStruct, "min")
##    mergedETStruct_max = merge_replicas(etStruct, "max")
##    print "merge replicas: min < median?", MA.add.reduce(MA.add.reduce(orng2ma(mergedETStruct_min[0][1][0]) > orng2ma(mergedETStruct_median[0][1][0]))) == 0
##    print "merge replicas: median < max?", MA.add.reduce(MA.add.reduce(orng2ma(mergedETStruct_median[0][1][0]) > orng2ma(mergedETStruct_max[0][1][0]))) == 0
##    print "merge replicas: abs(min - median)?", MA.add.reduce(MA.add.reduce(MA.absolute(orng2ma(mergedETStruct_mean[0][1][0]) - orng2ma(mergedETStruct_median[0][1][0]))))

