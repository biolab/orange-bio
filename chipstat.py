import math
import Numeric, MA
import LinearAlgebra
import stats
import NumExtn
import os, os.path
import orange

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
        merge_func = NumExtn.medianMA
    elif type == "min":
        merge_func = NumExtn.minMA
    elif type == "max":
        merge_func = NumExtn.maxMA
    else:
        raise AttributeError, "type = ['mean' | 'median' | 'min' | 'max']"
    for st, etList in aETStruct:
        shape[2] = len(etList)
        ma3d = MA.zeros(shape, MA.Float)
        for idx, et in enumerate(etList):
            ma3d[:,:,idx] = orng2ma(et)
        mergedETStruct.append((st, [ma2orng_keepClassMetas(merge_func(ma3d, 2), etList[0])]))
    return mergedETStruct
        


###########################################################################
#### average expression profile
###########################################################################
##
##def average_replicas(aETStruct):
##    """Returns a list of tuples (strain, avrg_orngET) where aETStruct corresponds to a list of tuples (strain, [orngET1, orngET2, ...]).
##    """
##    shape = [0,0,0]
##    et0 = aETStruct[0][1][0]                                            # the first example table
##    shape[0] = len(et0)                                                 # number of examples (genes)
##    shape[1] = len(et0.domain.attributes)                               # number of attributes (time points)
##    avrgETStruct = []
##    for st, etList in aETStruct:
##        shape[2] = len(etList)
##        ma3d = MA.zeros(shape, MA.Float)
##        for idx, et in enumerate(etList):
##            ma3d[:,:,idx] = orng2ma(et)
##        avrgETStruct.append((st, ma2orng(MA.average(ma3d, 2), et.domain)))
##    return avrgETStruct


#########################################################################
## ANOVA on genes
#########################################################################

class AnovaOnGeneResults:
    """A structure to store the results of ANOVA on genes where aETStruct corresponds to a list of tuples (strain, [orngET1, orngET2, ...]).
    """
    def __init__(self, aETStruct, anovaList):
        self.classNames = map(lambda x: x[0], aETStruct)
        self.anovaList = anovaList


def anova_on_genes(aETStruct, repMeasuresOnTime=False, callback = None):
    """Conducts ANOVA on genes and returns a tuple (list of [pTime, pStrain, pInter], list of ANOVA objects needed to conduct posthoc tests)
    where aETStruct corresponds to a list of tuples (strain, [orngET1, orngET2, ...]).
    Note: time points that cause empty cells are removed prior to conducting ANOVA.
    """
    ma3d = etStruct2ma3d(aETStruct)
    rgLens = Numeric.array(map(lambda x: len(x[1]), aETStruct))
    # arrays to store p-vals
    FprobTime = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
    FprobStrain = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
    FprobTimeStrain = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
    # store ANOVA objects
    anovaList = []
    # decide between non-repeated / repeated measures ANOVA for factor time
    if repMeasuresOnTime:
        fAnova = AnovaRM12LR
    else:
        fAnova = Anova2wayLR
    # check for empty cells for all genes at once and remove them
    tInd2rem = []
    ax2Ind = Numeric.concatenate(([0], Numeric.add.accumulate(rgLens)))
    for tIdx in range(ma3d.shape[1]):
        for rIdx in range(rgLens.shape[0]):
            if Numeric.add.reduce(MA.count(ma3d[:,tIdx,ax2Ind[rIdx]:ax2Ind[rIdx+1]],1)) == 0:
                tInd2rem.append(tIdx)
                break
    if len(tInd2rem) > 0:
        print "Warning: removing time indices %s for all genes" % (str(tInd2rem))
        tInd2keep = range(ma3d.shape[1])
        for tIdx in tInd2rem:
            tInd2keep.remove(tIdx)
        ma3d = MA.take(ma3d, tInd2keep, 1)
    # for each gene...
    for gIdx in range(ma3d.shape[0]):
        # check for empty cells for that gene
        tInd2rem = []
        for tIdx in range(ma3d.shape[1]):
            for rIdx in range(rgLens.shape[0]):
                if MA.count(ma3d[gIdx,tIdx,ax2Ind[rIdx]:ax2Ind[rIdx+1]]) == 0:
                    tInd2rem.append(tIdx)
                    break
        ma2d = ma3d[gIdx]
        if len(tInd2rem) > 0:
            print "Warning: removing time indices %s for gene %i" % (str(tInd2rem), gIdx)
            tInd2keep = range(ma3d.shape[1])
            for tIdx in tInd2rem:
                tInd2keep.remove(tIdx)
            ma2d = MA.take(ma2d, tInd2keep, 0)
        # conduct anova ANOVA
        an = fAnova(ma2d, rgLens, 1)
        if callback:
            callback()
        FprobTime[gIdx] = an.FAprob
        FprobStrain[gIdx] = an.FBprob
        FprobTimeStrain[gIdx] = an.FABprob
        anovaList.append(an)
    return Numeric.transpose(Numeric.array([FprobTime, FprobStrain, FprobTimeStrain])).tolist(), AnovaOnGeneResults(aETStruct, anovaList)


#########################################################################
## POSTHOC ANOVA on genes
#########################################################################

def posthoc_anova_on_genes(className1, className2, aAnovaOnGeneResults, callback=None):
    """Conduct posthoc anova test for the given class names and returns a list of of [pTime, pStrain, pInter].
    """
    cIdx1 = aAnovaOnGeneResults.classNames.index(className1)
    cIdx2 = aAnovaOnGeneResults.classNames.index(className2)
    # arrays to store p-vals
    numGenes = len(aAnovaOnGeneResults.anovaList)
    FpTime = -1*Numeric.ones((numGenes,), Numeric.Float)
    FpStrain = -1*Numeric.ones((numGenes,), Numeric.Float)
    FpTimeStrain = -1*Numeric.ones((numGenes,), Numeric.Float)
    for idx, an in enumerate(aAnovaOnGeneResults.anovaList):
        FpTime[idx], FpStrain[idx], FpTimeStrain[idx] = an.single_posthoc_anova_B_AB(cIdx1, cIdx2)
        if callback: callback()
    return Numeric.transpose(Numeric.array([FpTime, FpStrain, FpTimeStrain])).tolist()
        

###################################################################################
## HELPER FUNCTIONS
###################################################################################

def scaleStd(nm, axis=0):
    """Returns new masked numarray with values scaled by dividing by Median Absolute Difference: MAD = median(|val(i)-median(val(i))|)."""
    nm = NumExtn.swapaxesMA(MA.asarray(nm, MA.Float), 0, axis)
    return NumExtn.swapaxesMA(nm / NumExtn.stdMA(nm,0), 0, axis)

def scaleMad(nm, axis=0):
    """Returns new masked numarray with values scaled by dividing by Median Absolute Difference: MAD = median(|val(i)-median(val(i))|)."""
    nm = NumExtn.swapaxesMA(MA.asarray(nm, MA.Float), 0, axis)
    return NumExtn.swapaxesMA(nm / NumExtn.madMA(nm,0), 0, axis)


def centerMean(nm, axis=0):
    """Returns new masked numarray with values centered by subtracting Median."""
    nm = NumExtn.swapaxesMA(MA.asarray(nm, MA.Float), 0, axis)
    return NumExtn.swapaxesMA(nm - MA.average(nm,0), 0, axis)


def centerMed(nm, axis=0):
    """Returns new masked numarray with values centered by subtracting Median."""
    nm = NumExtn.swapaxesMA(MA.asarray(nm, MA.Float), 0, axis)
    return NumExtn.swapaxesMA(nm - NumExtn.medianMA(nm,0), 0, axis)


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
    ma3d = MA.zeros(shape, MA.Float)
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
    ma = MA.array(vals, MA.PyObject)
    if aExampleTable.domain.classVar != None:
        ma = ma[:,:-1]
    mask = MA.where(MA.equal(ma, "?"), 1, 0)
    for varIdx, var in enumerate(aExampleTable.domain.attributes):
        if type(var) != orange.FloatVariable:
            mask[:,varIdx] = Numeric.ones(len(aExampleTable))
    return MA.array(MA.array(ma, MA.PyObject, mask=mask).filled(1e20), MA.Float, mask=mask)


def ma2orng(arr2d, aDomain):
    """Converts MA.array to orange.ExampleTable where examples correspond to rows and attributes to columns.
    masked values converted to "?" (don't knows)
    arr2d.shape[1] must be equal to the number of the attributes of the given domain
    domain attributes mut be of type orange.FloatVariable
    """
    arr2d = MA.asarray(arr2d, MA.PyObject)
    assert MA.rank(arr2d) == 2, "2d array expected"
    assert len(aDomain.attributes) == arr2d.shape[1], "the shape of the array incompatible with the given domain"
    if aDomain.classVar != None:
        et = orange.ExampleTable(orange.Domain(aDomain.attributes, None), arr2d.tolist("?"))
        return orange.ExampleTable(aDomain, et)
    else:
        return orange.ExampleTable(aDomain, arr2d.tolist("?"))


def ma2orng_keepClassMetas(arr2d, aExampleTable):
    """Creates new example table where attribute values correspond to the given 2D array, class and meta attributes remain unchanged.
    """
    arr2d = MA.asarray(arr2d, MA.PyObject)
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
    return orange.ExampleTable([etAtt, etClassMeta])


###################################################################################
## PURE ANOVA BASED ON LINEAR REGRESSION
###################################################################################

class AnovaLRBase:
    """Base class for all ANOVA types based on multivariate linear regression: manipulation with dummy variables.
    """

    def getDummyRC(self, numVars, numReplicas=1):
        """Returns 2D array of reference cell encoded dummy variables.
        """
        if not type(numReplicas) == int:
            numReplicas = Numeric.asarray(numReplicas)
            assert len(numReplicas) == numVars, "numReplicas: int or argument of length numVars expected"
        return Numeric.repeat(Numeric.concatenate((Numeric.zeros((1,numVars-1), Numeric.Float), Numeric.identity(numVars-1, Numeric.Float)), 0), numReplicas, 0)
        

    def getDummyEE(self, numVars, numReplicas=1):
        """Returns 2D array of effects endoced dummy variables, shape (n, n-1) where n = numVars*numReplicas | sum(numReplicas).
        Each row represents a code for a single variable value, encoded by numVar-1 dummy variables.
        numReplicas = int | list | array of shape (numVars,): the number (or list) of observations per cell.
        """
        assert numVars > 0, "numVars <= 0"
        if type(numReplicas) == int:
            numReplicas = Numeric.ones((numVars,)) * numReplicas
        else:
            numReplicas = Numeric.asarray(numReplicas)
            assert len(numReplicas) == numVars, "numReplicas: int or argument of length numVars expected"
        if numVars == 1:
            return Numeric.zeros((numReplicas[0], 0))
        else:
            return Numeric.repeat(Numeric.concatenate((Numeric.identity(numVars-1, Numeric.Int), -1*Numeric.ones((1,numVars-1), Numeric.Int)), 0), numReplicas, 0)
        

    def getSubjectDummyEE(self, numSubjList, numReplicas=1):
        """Returns 2D array of effects endoced subject dummy variables, shape (sum(ni), sum(ni-1)) where numSubjList=[n0,n1,...].
        numSubjList: a list of the number of subjects within each level of the non-repeated factor.
        numReplicas: int | list | array of shape sum(ni): the number of times each subject code (row) is repeated.
        """
        numSubjList = Numeric.asarray(numSubjList)
        assert len(numSubjList.shape) == 1, "len(numSubjList.shape) != 1"
        if type(numReplicas) == int:
            numReplicas = Numeric.ones((numSubjList.shape[0],)) * numReplicas
        else:
            numReplicas = Numeric.asarray(numReplicas)
            assert len(numReplicas.shape) == 1, "len(numReplicas.shape) != 1"
            assert numReplicas.shape[0] == numSubjList.shape[0], "numReplicas: int or a list of the same length as numSubjList expected"
        si0 = Numeric.concatenate(([0], Numeric.add.accumulate(numSubjList)))      # indices for axis 0 of dm
        si1 = Numeric.concatenate(([0], Numeric.add.accumulate(numSubjList-1)))    # indices for axis 1 of dm
        dm = Numeric.zeros((si0[-1], si1[-1]))                                     # subject dummy codes matrix
        numRepeats = Numeric.ones((si0[-1],))
        for i in range(numSubjList.shape[0]):
            if numSubjList[i] == 1: print "Warning: a single subject in group %i" % (i)
            dm[si0[i]:si0[i+1], si1[i]:si1[i+1]] = self.getDummyEE(numSubjList[i],1)
            numRepeats[si0[i]:si0[i+1]] = Numeric.ones((numSubjList[i],)) * numReplicas[i]
        return Numeric.repeat(dm, numRepeats, 0)


    def getDummiesJoin(self, dummyMtrx1, dummyMtrx2):
        """Returns a tuple of resized dummy matrices containing the cartesian product of rows.
        """
        dummyMtrx1 = Numeric.repeat(dummyMtrx1, dummyMtrx2.shape[0], 0)                         # repeat each row separately
        dummyMtrx2 = Numeric.resize(dummyMtrx2, (dummyMtrx1.shape[0], dummyMtrx2.shape[1]))     # repeat whole matrix, add new rows
        return (dummyMtrx1, dummyMtrx2)


    def getDummyInteraction(self, dummyMtrx1, dummyMtrx2):
        """Returns a product of dummy matrices resized to the cartesian product of columns.
        """
        dmIntShp = (dummyMtrx1.shape[0], dummyMtrx1.shape[1] * dummyMtrx2.shape[1])
        dmInt1 = Numeric.repeat(dummyMtrx1, dummyMtrx2.shape[1], 1)                             # repeat each column separately
        dmInt2 = Numeric.transpose(Numeric.resize(Numeric.transpose(dummyMtrx2, (1,0)), (dmIntShp[1], dmIntShp[0])), (1,0)) # repeat whole matrix, add new columns
        return dmInt1*dmInt2


class Anova2wayLRBase(AnovaLRBase):
    """Base class for all 2 way ANOVA types (repeated and non-repeated measures) includes various post-hoc tests.
    The following variables should be defined by base classes:
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens
        self._MSpool = LR_treat.MSres
        self._DFpool = LR_treat.MSreg
    """

    def posthoc_tstat_B(self):
        """Fisher LSD: Compares all levels of factor B and returns (b,b) matrix of t-stat values and a matrix of corresponding p-values.
        Use when there is no significant interaction.
        Warning: Statistica computes it WRONG!
        """
        b = self._groupLens.shape[0]
        # precompute averages
        x_avrg = Numeric.zeros((b,), Numeric.Float)         # mean of cell means over factor A
        sum_count_x = Numeric.zeros((b,), Numeric.Float)
        groupLensAcc0 = Numeric.array([0] + Numeric.add.accumulate(self._groupLens).tolist())
        for idx in range(b):
            groupInd = range(groupLensAcc0[idx], groupLensAcc0[idx+1])
            xi = MA.take(self._arr2d, groupInd, 1)
            x_avrg[idx] = MA.average(MA.average(xi,1))                  # first average over replicas to obtain cell mean, then over factor A
            sum_count_x[idx] = Numeric.add.reduce(1./MA.count(xi,1))    # first count the number of measurements within cells, then sum inverses over factor A
            ##x_avrg[idx] = MA.average(MA.ravel(xi))                    # this is how Statistica computes it, but it is WRONG!
            ##sum_count_x[idx] = 1./MA.count(xi)                        # this is how Statistica computes it, but it is WRONG!
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (x_avrg.shape[0], x_avrg.shape[0]))
        sum_count_x_2 = Numeric.resize(sum_count_x, (sum_count_x.shape[0],sum_count_x.shape[0]))
        tstat = (MA.transpose(x_avrg_2) - x_avrg_2) * self._arr2d.shape[0] / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))
        ##tstat = (MA.transpose(x_avrg_2) - x_avrg_2) / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))     # this is how Statistica computes it, but it is WRONG!
        tprob = NumExtn.triangularPut(stats.abetai(0.5*self._DFpool, 0.5, float(self._DFpool) / (self._DFpool + NumExtn.triangularGet(tstat**2))),1,1)
        return tstat, tprob


    def isSignif_holm_B(self, alpha):
        """Conduct Holm step-down sequential multiple comparison on factor B.
        Returns (b,b) matrix indicating significant differences between levels of factor B and the direction of changes [-1|0|1].
        """
        tstat, pvals = self.posthoc_tstat_B()
        pvals = NumExtn.triangularGet(pvals)
        sortInd = Numeric.argsort(pvals)
        k = pvals.shape[0]
        isSignif = -1*Numeric.ones((k,), Numeric.Int)
        for j in range(k):
            isSignif[sortInd[j]] = pvals[sortInd[j]] < float(alpha) / (k-j)
        return NumExtn.triangularPut(isSignif, upper=1, lower=1) * (MA.greater(tstat,0) - MA.less(tstat,0))


    def posthoc_tstat_BbyA(self):
        """Fisher LSD: Compares all levels of factor B by each level of factor A separately.
        Returns: (b,b,a) matrix of t-stat values
                 (b,b,a) matrix of corresponding p-values
                 (b,b) matrix of geometric means of p-values over all levels of factor A.
        Use when there is a significant interaction and factor A is of no interest.
        """
        a = self._arr2d.shape[0]
        b = self._groupLens.shape[0]
        # precompute averages
        x_avrg = Numeric.zeros((a,b), Numeric.Float)
        sum_count_x = Numeric.zeros((a,b), Numeric.Float)
        groupLensAcc0 = Numeric.array([0] + Numeric.add.accumulate(self._groupLens).tolist())
        for idx in range(b):
            groupInd = range(groupLensAcc0[idx], groupLensAcc0[idx+1])
            xi = MA.take(self._arr2d, groupInd, 1)
            x_avrg[:,idx] = MA.average(xi, 1)
            sum_count_x[:,idx] = 1. / MA.count(xi, 1)
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (b,a,b))
        sum_count_x_2 = Numeric.resize(sum_count_x, (b,a,b))
        tstat = (MA.transpose(x_avrg_2, (2,1,0)) - x_avrg_2) / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2, (2,1,0)) + sum_count_x_2))
        # get p-values for each level of factor a separately (axis 1)
        tprob = MA.array(-1*Numeric.ones(tstat.shape, Numeric.Float), mask=Numeric.ones(tstat.shape))
        for idx1 in range(tprob.shape[1]):
            tprob[:,idx1,:] = NumExtn.triangularPut(stats.abetai(0.5*self._DFpool, 0.5, float(self._DFpool) / (self._DFpool + NumExtn.triangularGet(tstat[:,idx1,:]**2))),1,1)
        return tstat, tprob


    def posthoc_anova_B_AB(self):
        """Conduct ANOVA for each pair of levels of factor B to detect interaction effect.
        Returns for each pair of factor B levels:
            (b,b) matrix of F-stat values for effect of factor B 
            (b,b) matrix of p-values for effect of factor B
            (b,b) matrix of F-stat values for interaction effect
            (b,b) matrix of p-values for interaction effect
        """
        if self._addInteraction != 1:
            raise "Error: posthoc_anova_B_AB can be conducted only when the interaction effect has been tested"
        b = self._groupLens.shape[0]
        FB = MA.masked * MA.ones((b,b), MA.Float)
        FBprob = MA.masked * MA.ones((b,b), MA.Float)
        FAB = MA.masked * MA.ones((b,b), MA.Float)
        FABprob = MA.masked * MA.ones((b,b), MA.Float)
        groupLensAcc0 = Numeric.array([0] + Numeric.add.accumulate(self._groupLens).tolist())
        groupInd = map(lambda i,j: range(i,j), groupLensAcc0[:-1],groupLensAcc0[1:])
        for i in range(b):
            for j in range(i+1,b):
                takeInd = groupInd[i] + groupInd[j]
                an = self.__class__(MA.take(self._arr2d, takeInd, 1), [self._groupLens[i],self._groupLens[j]], addInteraction=1)
                FB[i,j] = FB[j,i] = an.FB
                FBprob[i,j] = FBprob[j,i] = an.FBprob
                FAB[i,j] = FAB[j,i] = an.FAB
                FABprob[i,j] = FABprob[j,i] = an.FABprob
        return FB, FBprob, FAB, FABprob


    def single_posthoc_anova_B_AB(self, bLevel1, bLevel2):
        """Conduct ANOVA for a specific pair of factor B levels and returns (FprobA, FprobB, FprobAB)
        """
        if self._addInteraction != 1:
            raise "Error: single_posthoc_anova_B_AB can be conducted only when the interaction effect has been tested"
        groupLensAcc0 = Numeric.array([0] + Numeric.add.accumulate(self._groupLens).tolist())
        groupInd = map(lambda i,j: range(i,j), groupLensAcc0[:-1],groupLensAcc0[1:])
        takeInd = groupInd[bLevel1] + groupInd[bLevel2]
        an = self.__class__(MA.take(self._arr2d, takeInd, 1), [self._groupLens[bLevel1],self._groupLens[bLevel2]], addInteraction=1)
        return an.FAprob, an.FBprob, an.FABprob

    
class Anova2wayLR(Anova2wayLRBase):
    """2 way ANOVA with multiple observations per cell for unbalanced designs using multivariate linear regression.
    Multiple observations given at axis 1 together with the list of group lenghts, e.g. [3,2].
    Supports balanced and unbalanced designs, i.e. different number of observations per cell and missing (masked) data.
    TODO: account for empty cells (this implementation results in singular design matrix)
    """

    def __init__(self, arr2d, groupLens, addInteraction=0):
        """arr2d: 2D masked array where replicas are given at axis 1 together with lenghts of the groups (multiple observations per cell).
        addInteraction: include / exclude the interaction between factors
        """
        arr2d = MA.asarray(arr2d)
        groupLens = Numeric.array(groupLens)    # make a copy for the case if subjects or levels of factor B are removed
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"

        # check if there exist factor-levels with all values missing (for A and B)
        # if so, reduce the corresponding DF and fix the number of dummy variables
        missIndA = Numeric.compress(MA.count(arr2d, 1) == 0, Numeric.arange(arr2d.shape[0]))
        if missIndA.shape[0] > 0:
            print "Warning: removig factor A level(s) %s" % str(missIndA.tolist())
            takeIndA = Numeric.compress(MA.count(arr2d, 1) != 0, Numeric.arange(arr2d.shape[0]))
            arr2d = MA.take(arr2d, takeIndA, 0)
        missIndSubj = Numeric.compress(MA.count(arr2d, 0) == 0, Numeric.arange(arr2d.shape[1]))
        if missIndSubj.shape[0] > 0:
            takeIndSubj = Numeric.compress(MA.count(arr2d, 0) != 0, Numeric.arange(arr2d.shape[1]))
            arr2d = MA.take(arr2d, takeIndSubj, 1)
            # fix groupLens
            mapSubj2Group = -1*Numeric.ones(arr2d.shape[1])
            groupLensAcc0 = [0] + Numeric.add.accumulate(groupLens).tolist()
            for i in range(len(groupLensAcc0) - 1):
                mapSubj2Group[groupLensAcc0[i]:groupLensAcc0[i+1]] = i
            for subjIdx in missIndSubj:
                groupLens[mapSubj2Group[subjIdx]] -= 1
        # fix number of factor B levels
        missIndB = Numeric.compress(groupLens <= 0, Numeric.arange(groupLens.shape[0]))
        if missIndB.shape[0] > 0:
            print "Warning: removig factor B level(s) %s" % str(missIndB.tolist())
            takeIndB = Numeric.compress(groupLens > 0, Numeric.arange(groupLens.shape[0]))
            groupLens = Numeric.take(groupLens, takeIndB)
            # arr2d has already been taken care of by removing the subjects without observations

        # check if data is balanced
        isBalanced = (1. * Numeric.add.reduce((1.*groupLens/groupLens[0]) == Numeric.ones(groupLens.shape[0])) / groupLens.shape[0]) == 1

        # check for empty cells, raise exception if empty cells exist and addInteraction=1
        if addInteraction:
            ax1Ind = Numeric.concatenate(([0], Numeric.add.accumulate(groupLens)))
            for idx in range(groupLens.shape[0]):
                if Numeric.add.reduce(MA.count(arr2d[:,ax1Ind[idx]:ax1Ind[idx+1]],1) == 0) > 0:
                    raise ValueError, "arr2d has empty cells, cannot evaluate the interaction between factors, try addInteraction=0"

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d) # [a0b0, a0b1, .., a1b0, ...]
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = NumExtn.compressIndices(self.y)

        # dummy variables
        self.dummyA = self.getDummyEE(arr2d.shape[0], 1)
        self.dummyB = self.getDummyEE(len(groupLens), groupLens)
        self.dummyA, self.dummyB = self.getDummiesJoin(self.dummyA, self.dummyB)
        if addInteraction:
            self.dummyAB = self.getDummyInteraction(self.dummyA, self.dummyB)
            self.dummyAB = Numeric.take(self.dummyAB, takeInd, 0)
        self.dummyA = Numeric.take(self.dummyA, takeInd, 0)
        self.dummyB = Numeric.take(self.dummyB, takeInd, 0)

        # run analysis
        if addInteraction:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB, self.dummyAB), 1), self.y)                
        else:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB),1), self.y)
            
        if isBalanced and noMissing:
            LR_A = MultLinReg(self.dummyA, self.y)
            LR_B = MultLinReg(self.dummyB, self.y)
            self.FA = LR_A.MSreg / LR_treat.MSres
            self.FB = LR_B.MSreg / LR_treat.MSres
            if addInteraction:
                LR_AB = MultLinReg(self.dummyAB, self.y)
                self.FAB = LR_AB.MSreg / LR_treat.MSres
        else:
            if addInteraction:
                LR_A_B = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB), 1), self.y)
                LR_A_AB = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyAB), 1), self.y)
                LR_B_AB = MultLinReg(Numeric.concatenate((self.dummyB, self.dummyAB), 1), self.y)
                self.FA = (LR_treat.SSreg - LR_B_AB.SSreg) / self.dummyA.shape[1] / LR_treat.MSres
                self.FB = (LR_treat.SSreg - LR_A_AB.SSreg) / self.dummyB.shape[1] / LR_treat.MSres
                self.FAB = (LR_treat.SSreg - LR_A_B.SSreg) / self.dummyAB.shape[1] / LR_treat.MSres
            else:
                LR_A = MultLinReg(self.dummyA, self.y)
                LR_B = MultLinReg(self.dummyB, self.y)
                self.FA = (LR_treat.SSreg - LR_B.SSreg) / self.dummyA.shape[1] / LR_treat.MSres
                self.FB = (LR_treat.SSreg - LR_A.SSreg) / self.dummyB.shape[1] / LR_treat.MSres

        self.FAprob = stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        self.FBprob = stats.fprob(self.dummyB.shape[1], LR_treat.DFres, self.FB)
        if addInteraction:
            self.FABprob = stats.fprob(self.dummyAB.shape[1], LR_treat.DFres, self.FAB)

        # store variables needed for posthoc tests
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens
        self._MSpool = LR_treat.MSres
        self._DFpool = LR_treat.MSreg


class AnovaRM12LR(Anova2wayLRBase):
    """2 way ANOVA with REPEATED MEASURES on factor A using multivariate linear regression.
    Axis 0 correspond to different levels of factor A, subjects given at axis 1, subjects nested inside factor B according to the given group lenghts.
    Factor A is a within-subject effect, factor B is a between-ubject effect.
    Example:
        arr2d = [[a0b0 a0b0 a0b1 a0b1 a0b1], [a1b0 a1b0 a1b1 a1b1 a1b1]]: 2 levels of factor A, 2 levels of factor B, 5 subjects
        groupLens = [2,3]
    Supports balanced and unbalanced designs, i.e. different number of observations per cell and missing (masked) data.
    TODO: fix the F and DF for factor B when there are missing values (see Glantz, pp. 462-5)
    """

    def __init__(self, arr2d, groupLens, addInteraction=0):
        """arr2d: 2D masked array: [[a0b0 a0b0 a0b1 a0b1 a0b1], [a1b0 a1b0 a1b1 a1b1 a1b1]]
        groupLens: lenghts of axis 1 indices that belong to the same level of factor B, e.g. [2,3]
        addInteraction: include / exclude the interaction between factors A and B
        """
        arr2d = MA.asarray(arr2d)
        groupLens = Numeric.array(groupLens)    # make a copy for the case if subjects or levels of factor B are removed
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"
        assert arr2d.shape[1] == Numeric.add.reduce(groupLens), "arr2d.shape[1] != Numeric.add.reduce(groupLens)"

        # check if there exist factor-levels with all values missing (for A, B and Subj)
        # if so, reduce the corresponding DF and fix the number of dummy variables
        missIndA = Numeric.compress(MA.count(arr2d, 1) == 0, Numeric.arange(arr2d.shape[0]))
        if missIndA.shape[0] > 0:
            print "Warning: removig factor A level(s) %s" % str(missIndA.tolist())
            takeIndA = Numeric.compress(MA.count(arr2d, 1) != 0, Numeric.arange(arr2d.shape[0]))
            arr2d = MA.take(arr2d, takeIndA, 0)
        missIndSubj = Numeric.compress(MA.count(arr2d, 0) == 0, Numeric.arange(arr2d.shape[1]))
        if missIndSubj.shape[0] > 0:
            print "Warning: removig subject(s) %s" % str(missIndSubj.tolist())
            takeIndSubj = Numeric.compress(MA.count(arr2d, 0) != 0, Numeric.arange(arr2d.shape[1]))
            arr2d = MA.take(arr2d, takeIndSubj, 1)
            # fix groupLens
            mapSubj2Group = -1*Numeric.ones(arr2d.shape[1])
            groupLensAcc0 = [0] + Numeric.add.accumulate(groupLens).tolist()
            for i in range(len(groupLensAcc0) - 1):
                mapSubj2Group[groupLensAcc0[i]:groupLensAcc0[i+1]] = i
            for subjIdx in missIndSubj:
                groupLens[mapSubj2Group[subjIdx]] -= 1
        # fix number of factor B levels
        missIndB = Numeric.compress(groupLens <= 0, Numeric.arange(groupLens.shape[0]))
        if missIndB.shape[0] > 0:
            print "Warning: removig factor B level(s) %s" % str(missIndB.tolist())
            takeIndB = Numeric.compress(groupLens > 0, Numeric.arange(groupLens.shape[0]))
            groupLens = Numeric.take(groupLens, takeIndB)
            # arr2d has already been taken care of by removing the subjects without observations

        # remove other missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d) # [a0b0, a0b1, .., a1b0, ...]
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = NumExtn.compressIndices(self.y)

        # check if data is balanced
        isBalanced = (1. * Numeric.add.reduce((1.*groupLens/groupLens[0]) == Numeric.ones(groupLens.shape[0])) / groupLens.shape[0]) == 1

        # dummy variables
        self.dummyA = self.getDummyEE(arr2d.shape[0], 1)
        self.dummyB = self.getDummyEE(groupLens.shape[0], groupLens)
        self.dummySubjInB = self.getSubjectDummyEE(groupLens, 1)
        dummyB_Subj = Numeric.concatenate((self.dummyB, self.dummySubjInB), 1)
        self.dummyA, dummyB_Subj = self.getDummiesJoin(self.dummyA, dummyB_Subj)
        self.dummyB = dummyB_Subj[:, 0:self.dummyB.shape[1]]
        self.dummySubjInB = dummyB_Subj[:, self.dummyB.shape[1]:]
        if addInteraction:
            self.dummyAB = self.getDummyInteraction(self.dummyA, self.dummyB)
            self.dummyAB = Numeric.take(self.dummyAB, takeInd, 0)
        self.dummyA = Numeric.take(self.dummyA, takeInd, 0)
        self.dummyB = Numeric.take(self.dummyB, takeInd, 0)
        self.dummySubjInB = Numeric.take(self.dummySubjInB, takeInd, 0)

        # run analysis
        if addInteraction:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA, self.dummyB, self.dummyAB),1), self.y)
        else:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA, self.dummyB),1), self.y)            

        if isBalanced and noMissing and False:
            # non-repeated-measures factor (B)
            LR_B = MultLinReg(self.dummyB, self.y)
            LR_SubjInB = MultLinReg(self.dummySubjInB, self.y)
            self.FB = LR_B.MSreg / LR_SubjInB.MSreg
            # repeated measures factor (A)
            LR_A = MultLinReg(self.dummyA, self.y)
            self.FA = LR_A.MSreg / LR_treat.MSres
            # interaction (AB)
            if addInteraction:
                LR_AB = MultLinReg(self.dummyAB, self.y)
                self.FAB = LR_AB.MSreg / LR_treat.MSres
            # store variables needed for posthoc tests
            MS_SubjInB = LR_SubjInB.MSreg
        else:
            ######################################################################################################################################
            ## print "Warning: missing data with 2 way repeated measures ANOVA: F and DF should be adjusted due to random effect factor (subjects)"
            ######################################################################################################################################
            if addInteraction:
                # non-repeated-measures factor (B)
                LR_S_A_AB = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA, self.dummyAB),1), self.y)
                LR_A_B_AB = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB, self.dummyAB),1), self.y)
                self.FB = (LR_treat.SSreg - LR_S_A_AB.SSreg) / self.dummyB.shape[1] / (LR_treat.SSreg - LR_A_B_AB.SSreg) * self.dummySubjInB.shape[1]
                # repeated measures factor (A)
                LR_S_B_AB = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyB, self.dummyAB),1), self.y)
                self.FA = (LR_treat.SSreg - LR_S_B_AB.SSreg) / self.dummyA.shape[1] / LR_treat.MSres
                # interaction (AB)
                LR_S_A_B = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA, self.dummyB),1), self.y)
                self.FAB = (LR_treat.SSreg - LR_S_A_B.SSreg) / self.dummyAB.shape[1] / LR_treat.MSres
                # store variables needed for posthoc tests
                MS_SubjInB = (LR_treat.SSreg - LR_A_B_AB.SSreg) / self.dummySubjInB.shape[1]
            else:
                # non-repeated-measures factor (B)
                LR_S_A = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA),1), self.y)
                LR_A_B = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB),1), self.y)
                self.FB = (LR_treat.SSreg - LR_S_A.SSreg) / self.dummyB.shape[1] / (LR_treat.SSreg - LR_A_B.SSreg) * self.dummySubjInB.shape[1]
                # repeated measures factor (A)
                LR_S_B = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyB),1), self.y)
                self.FA = (LR_treat.SSreg - LR_S_B.SSreg) / self.dummyA.shape[1] / LR_treat.MSres
                # store variables needed for posthoc tests
                MS_SubjInB = (LR_treat.SSreg - LR_A_B.SSreg) / self.dummySubjInB.shape[1]

        self.FBprob = stats.fprob(self.dummyB.shape[1], self.dummySubjInB.shape[1], self.FB)
        self.FAprob = stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        if addInteraction:
            self.FABprob = stats.fprob(self.dummyAB.shape[1], LR_treat.DFres, self.FAB)

        # store variables needed for posthoc tests
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens
        # estimate the variance and degrees of freedom by pooling residual variance for factors A and B
        MS_res = LR_treat.MSres
        a = self.dummyA.shape[1] + 1    # number of levels of factor A
        DF_SubjInB = self.dummySubjInB.shape[1] # = LR_SubjInB.DFreg
        DF_res = LR_treat.DFres
        self._MSpool = (MS_SubjInB + (a-1)*MS_res) / a
        self._DFpool = ((MS_SubjInB + (a-1)*MS_res)**2) / ((MS_SubjInB**2 / DF_SubjInB) + (((a-1)*MS_res)**2 / DF_res))


class MultLinReg:
    """Multivariate Linear Regression.
    """
    def __init__(self, X, y, summary=0):
        """
        x[i,j]: ith observation for the jth independent variable, shape (n,k)
        y[i]: ith observations of the dependent variable, shape (n,)
        """
        X = Numeric.asarray(X)
        y = Numeric.asarray(y)
        assert len(X.shape) == 2, "len(X.shape) != 2"
        assert len(y.shape) == 1, "len(y.shape) != 1"
        self.DFreg = X.shape[1]                         # number of independent variables (k)
        self.DFres = X.shape[0] - (X.shape[1]) - 1      # number of observations (n) - number of ind. variables (k) - 1
        self.DFtot = X.shape[0] - 1                     # number of observations (n) - 1
        X = Numeric.concatenate((Numeric.ones((X.shape[0],1), Numeric.Float), X), 1)
        XT = Numeric.transpose(X)
        try:
            self._XTXinv = LinearAlgebra.inverse(Numeric.dot(XT,X))
        except LinearAlgebra.LinAlgError:
            print "Warning: singular matrix, using generalized_inverse"
            self._XTX = Numeric.dot(XT,X)   # store the singuar matrix
            self._XTXinv = LinearAlgebra.generalized_inverse(self._XTX)
        y_mean = Numeric.average(y)
        # regression
        self.b = Numeric.dot(Numeric.dot(self._XTXinv, XT), y)
        self.y_hat = Numeric.dot(X, self.b)
        # release the variables not needed
        X = None
        XT = None
        # SS
        self.SSreg = Numeric.add.reduce((self.y_hat - y_mean)**2)
        self.MSreg = self.SSreg / self.DFreg
        self.SSres = Numeric.add.reduce((y - self.y_hat)**2)
        self.MSres = self.SSres / self.DFres                    # equals to square std. error: s^2_{y|x}
        if self.MSres == 0.0:
            print "Warning MultLinReg: MSres equals 0, replaced by 1e-20"
            self.MSres = 1e-20
        self.SStot = self.SSreg + self.SSres                    # equal to: self.SStot = Numeric.add.reduce((y - y_mean)**2)
        self.MStot = self.SStot / self.DFtot
        # regression summary
        if summary: self.getRegressionSummary()


    def __call__(self, X):
        X = Numeric.asarray(X)
        if len(X.shape) == 1:
            X = Numeric.concatenate(([1], X),0)
        elif len(X.shape) == 2:
            X = Numeric.concatenate((Numeric.ones((X.shape[0],1), Numeric.Float), X), 1)
        return Numeric.dot(X, self.b)

 
    def getRegressionSummary(self):
        self.stdErr = math.sqrt(self.MSres)                         # std. error of the regression estimate: s_{x|y}
        self.VarCovar = self.MSres * self._XTXinv                   # variance-covariance matrix: s^2_b
        self.b_se = Numeric.sqrt(Numeric.diagonal(self.VarCovar))   # std. error of regression coef.: s_b_i
        self.b_t = self.b / self.b_se                               # t-test whether each regression coef. is significantly different from 0
        self.b_tprob = stats.abetai(0.5*self.DFres, 0.5, float(self.DFres)/(self.DFres + self.b_t**2))
        # correlation between coefficients (b's)
        # ERROR in Glantz et al.: self.Corr = Numeric.sqrt(self.VarCovar) / Numeric.sqrt(Numeric.dot(self.b_se[:,Numeric.NewAxis],self.b_se[Numeric.NewAxis,:]))
        self.Corr = self.VarCovar / Numeric.dot(self.b_se[:,Numeric.NewAxis],self.b_se[Numeric.NewAxis,:])
        self.F = self.MSreg / self.MSres        # overall goodness of fit 
        self.Fprob = stats.fprob(self.DFreg, self.DFres, self.F)
        self.R2 = self.SSreg/self.SStot         # coef. of determination: fraction of variance explained by regression plane
        self.R2_adj = 1-(self.MSres/self.MStot) # adjusted coef. of determination (unbiased estimator of population R2; larget the better)
        self.R = math.sqrt(self.R2)             # multiple correlation coefficeint: how tight the data points cluster around the regression plane



if __name__ == "__main__":
    
    import Meda.Preproc
    reload(Meda.Preproc)
    import Dicty
    reload(Dicty)
    import Dicty.DAnnotation

##    myPath = r"C:\Python23\Lib\site-packages\orange\OrangeWidgets\Genomics\chipdata"
##    myPath = r"C:\Documents and Settings\peterjuv\My Documents\Transcriptional Phenotype\DictyData\orngDataStruct_Nancy"
    myPath = r"C:\Documents and Settings\peterjuv\My Documents\Transcriptional Phenotype\DictyData\orngDataStruct_Nancy_yakApufA"

    def generateETStruct(path):
        ddbList = Dicty.DAnnotation.getDDBList()
        if not os.path.exists(path):
            os.mkdir(path)
        DN = Dicty.DData.DData_Nancy()    
        for st in DN.strains:
            pathSt = path + "\\" + st
            if not os.path.exists(pathSt):
                os.mkdir(pathSt)
            for rep in DN.strain2replicaList(st):
                ma2d = DN.getRaw2d(rep)
                et = Meda.Preproc.ma2orng(ma2d,Meda.Preproc.getTcDomain(ma2d.shape[1], False, [], None))
                et.domain.addmeta(orange.newmetaid(), orange.StringVariable("DDB"))
                for eIdx,e in enumerate(et):
                    e["DDB"] = ddbList[eIdx]
                orange.saveTabDelimited(pathSt + "\\" + rep + ".tab", et)

    # generate orange example table structure
##    generateETStruct(myPath)

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

    # test Anova
##    pvals, aAnovaOnGeneResults = anova_on_genes(etStruct)
##    anova_meda = Dicty.DAnova.DGeneAnovaSTreduced(ma3d, DNpufAyakA.replicaGroupInd, 0.05, 0)
##    pvals_meda = Numeric.transpose(Numeric.array([myAnova.FprobTime, myAnova.FprobStrain, myAnova.FprobTimeStrain]))
##    pvals_OW = Numeric.array(pvals)
##    print "Anova OK?", Numeric.add.reduce(Numeric.add.reduce(pvals_meda - pvals_OW)) < 0.00001
    
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

