## Automatically adapted for numpy.oldnumeric Oct 03, 2007 by 

from __future__ import absolute_import

import math

import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA
import numpy.oldnumeric.linear_algebra as LinearAlgebra

import scipy.stats

from . import numpyExtn

#######################################################################################
## ANOVA based on Multivariate Linear Regression
#######################################################################################

class AnovaLRBase:
    """Base class for all ANOVA types based on multivariate linear regression: manipulation with dummy variables.
    """

    def getDummyRC(self, numVars, numReplicas=1):
        """Returns 2D array of reference cell encoded dummy variables.
        """
        if not isinstance(numReplicas, int):
            numReplicas = Numeric.asarray(numReplicas)
            assert len(numReplicas) == numVars, "numReplicas: int or argument of length numVars expected"
        return Numeric.repeat(Numeric.concatenate((Numeric.zeros((1,numVars-1), Numeric.Float), Numeric.identity(numVars-1, Numeric.Float)), 0), numReplicas, 0)
        

    def getDummyEE(self, numVars, numReplicas=1):
        """Returns 2D array of effects endoced dummy variables, shape (n, n-1) where n = numVars*numReplicas | sum(numReplicas).
        Each row represents a code for a single variable value, encoded by numVar-1 dummy variables.
        numReplicas = int | list | array of shape (numVars,): the number (or list) of observations per cell.
        """
        assert numVars > 0, "numVars <= 0"
        if isinstance(numReplicas, int):
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
        if isinstance(numReplicas, int):
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
            xi = self._arr2d.take(groupInd, 1)
            x_avrg[idx] = MA.average(MA.average(xi,1))                  # first average over replicas to obtain cell mean, then over factor A
            sum_count_x[idx] = Numeric.add.reduce(1./MA.count(xi,1))    # first count the number of measurements within cells, then sum inverses over factor A
            ##x_avrg[idx] = MA.average(MA.ravel(xi))                    # this is how Statistica computes it, but it is WRONG!
            ##sum_count_x[idx] = 1./MA.count(xi)                        # this is how Statistica computes it, but it is WRONG!
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (x_avrg.shape[0], x_avrg.shape[0]))
        sum_count_x_2 = Numeric.resize(sum_count_x, (sum_count_x.shape[0],sum_count_x.shape[0]))
        tstat = (MA.transpose(x_avrg_2) - x_avrg_2) * self._arr2d.shape[0] / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))
        ##tstat = (MA.transpose(x_avrg_2) - x_avrg_2) / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))     # this is how Statistica computes it, but it is WRONG!
        tprob = numpyExtn.triangularPut(scipy.stats.abetai(0.5*self._DFpool, 0.5, float(self._DFpool) / (self._DFpool + numpyExtn.triangularGet(tstat**2))),1,1)
        return tstat, tprob


    def isSignif_holm_B(self, alpha):
        """Conduct Holm step-down sequential multiple comparison on factor B.
        Returns (b,b) matrix indicating significant differences between levels of factor B and the direction of changes [-1|0|1].
        """
        tstat, pvals = self.posthoc_tstat_B()
        pvals = numpyExtn.triangularGet(pvals)
        sortInd = Numeric.argsort(pvals)
        k = pvals.shape[0]
        isSignif = -1*Numeric.ones((k,), Numeric.Int)
        for j in range(k):
            isSignif[sortInd[j]] = pvals[sortInd[j]] < float(alpha) / (k-j)
        return numpyExtn.triangularPut(isSignif, upper=1, lower=1) * (MA.greater(tstat,0) - MA.less(tstat,0))


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
            xi = self._arr2d.take(groupInd, 1)
            x_avrg[:,idx] = MA.average(xi, 1)
            sum_count_x[:,idx] = 1. / MA.count(xi, 1)
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (b,a,b))
        sum_count_x_2 = Numeric.resize(sum_count_x, (b,a,b))
        tstat = (MA.transpose(x_avrg_2, (2,1,0)) - x_avrg_2) / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2, (2,1,0)) + sum_count_x_2))
        # get p-values for each level of factor a separately (axis 1)
        tprob = MA.array(-1*Numeric.ones(tstat.shape, Numeric.Float), mask=Numeric.ones(tstat.shape))
        for idx1 in range(tprob.shape[1]):
            tprob[:,idx1,:] = numpyExtn.triangularPut(scipy.stats.abetai(0.5*self._DFpool, 0.5, float(self._DFpool) / (self._DFpool + numpyExtn.triangularGet(tstat[:,idx1,:]**2))),1,1)
        return tstat, tprob


    def posthoc_anova(self):
        """Conduct ANOVA for each pair of levels of factor B.
        Returns for each pair of factor B levels:
            (b,b) matrix of F-stat values for effect of factor A
            (b,b) matrix of p-values for effect of factor A
            (b,b) matrix of F-stat values for effect of factor B 
            (b,b) matrix of p-values for effect of factor B
            (b,b) matrix of F-stat values for interaction effect
            (b,b) matrix of p-values for interaction effect
        """
        if self._addInteraction != 1:
            raise "Error: posthoc_anova can be conducted only when the interaction effect has been tested"
        b = self._groupLens.shape[0]
        FA = MA.masked * MA.ones((b,b), Numeric.Float)
        FAprob = MA.masked * MA.ones((b,b), Numeric.Float)
        FB = MA.masked * MA.ones((b,b), Numeric.Float)
        FBprob = MA.masked * MA.ones((b,b), Numeric.Float)
        FAB = MA.masked * MA.ones((b,b), Numeric.Float)
        FABprob = MA.masked * MA.ones((b,b), Numeric.Float)
        groupLensAcc0 = Numeric.array([0] + Numeric.add.accumulate(self._groupLens).tolist())
        groupInd = map(lambda i,j: range(i,j), groupLensAcc0[:-1],groupLensAcc0[1:])
        for i in range(b):
            for j in range(i+1,b):
                takeInd = groupInd[i] + groupInd[j]
                an = self.__class__(self._arr2d.take(takeInd, 1), [self._groupLens[i],self._groupLens[j]], addInteraction=1)
                FA[i,j] = FA[j,i] = an.FA
                FAprob[i,j] = FAprob[j,i] = an.FAprob
                FB[i,j] = FB[j,i] = an.FB
                FBprob[i,j] = FBprob[j,i] = an.FBprob
                FAB[i,j] = FAB[j,i] = an.FAB
                FABprob[i,j] = FABprob[j,i] = an.FABprob
        return FA, FAprob, FB, FBprob, FAB, FABprob


class Anova1wayLR(AnovaLRBase):
    """1 way ANOVA for unbalanced designs using multivariate linear regression.
    Multiple observations given as a list of indices groups, e.g. [[0,1,2],[3,4]].
    BUG: the indices are not taken into account, but only the lenghts of the groups!
    Supports balanced and unbalanced designs and missing (masked) data.
    """

    def __init__(self, arr1d, replicaGroupLens):
        """arr1d: 1D masked array: [x1,x2,x3, y1,y2, z1,z2,z3,z4] where x,y,z correspond to factor-levels with 3,2 and 4 replicas, respectively
        replicaGroupLens: the number of replicas in individual groups, i.e. [3,2,4]
        """
        arr1d = MA.asarray(arr1d)
        assert len(arr1d.shape) == 1, "len(arr1d.shape) != 1"

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.asarray(arr1d)
        self.y, takeInd = numpyExtn.compressIndices(self.y)

        # dummy variables
        self.dummyA = self.getDummyEE(len(replicaGroupLens), replicaGroupLens)
        self.dummyA = Numeric.take(self.dummyA, takeInd, 0)

        # run analysis
        LRA = MultLinReg(self.dummyA, self.y)
        try:
            self.F = LRA.MSreg / LRA.MSres
            self.Fprob = scipy.stats.fprob(LRA.DFreg, LRA.DFres, self.F)
        except ZeroDivisionError:
            self.F = 0
            self.Fprob = 1


class Anova1wayLR_2D(AnovaLRBase):
    """1 way ANOVA for unbalanced designs using multivariate linear regression.
    Axis 0 correspond to replicas and axis 1 correspond to different factor-levels.
    Supports balanced and unbalanced designs and missing (masked) data.
    """

    def __init__(self, arr2d):
        arr2d = MA.asarray(arr2d)
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(MA.transpose(arr2d))
        self.y, takeInd = numpyExtn.compressIndices(self.y)

        # adjust degrees of freedom for factor-leves that contain no data
        numFactorLevels = int(MA.add.reduce(MA.greater(arr2d.count(0), 0)))
        zeroLengthFactorLevels = Numeric.compress(MA.equal(arr2d.count(0), 0), range(arr2d.shape[1]))
        if numFactorLevels > 1:

##            # dummy variables (BUG: index error when the number of levels is reduced)
##            print numFactorLevels, arr2d.shape[0], takeInd
##            self.dummyA = self.getDummyEE(numFactorLevels, arr2d.shape[0])
##            self.dummyA = Numeric.take(self.dummyA, takeInd, 0)
            # dummy variables
            self.dummyA = self.getDummyEE(arr2d.shape[1], arr2d.shape[0])   # create all dummies
            for lev in zeroLengthFactorLevels:                              # reduce length of dummies to adjust for actual number of factor A levels
                self.dummyA = Numeric.concatenate([self.dummyA[:,0:lev],self.dummyA[:,lev+1:]],1)
            if len(zeroLengthFactorLevels) > 0 and zeroLengthFactorLevels[-1] == arr2d.shape[1]-1:
                self.dummyA[self.dummyA.shape[0]-arr2d.shape[1]-1:] = -1    # if the last factor level is missing, manually set the last dummies to matrix of -1's
            self.dummyA = Numeric.take(self.dummyA, takeInd, 0)             # take dummies where data is not masked

            # run analysis
            LRA = MultLinReg(self.dummyA, self.y)
            try:
                self.F = LRA.MSreg / LRA.MSres
                self.Fprob = scipy.stats.fprob(LRA.DFreg, LRA.DFres, self.F)
            except ZeroDivisionError:
                self.F = 0
                self.Fprob = 1

        else:
            print "Giving up ANOVA: a factor with a single level encountered."
            self.F = 0
            self.Fprob = 1


class _Anova2wayLR_bug_missing_factor_leves(AnovaLRBase):
    """2 way ANOVA with multiple observations per cell for unbalanced designs using multivariate linear regression.
    Multiple observations given at axis 1 together with the list of indices groups, e.g. [[0,1,2],[3,4]].
    Supports balanced and unbalanced designs, i.e. different number of observations per cell and missing (masked) data.
    TODO: account for empty cells (this implementation results in singular design matrix)
    """

    def __init__(self, arr2d, replicaGroupInd, addInteraction=0):
        """arr2d: 2D masked array where replicas are given at axis 1 together with the indices groups.
        replicaGroupInd: indices of axis 1 that belong to a common factor (multiple observations per cell)
        addInteraction: include / exclude the interaction between factors
        """
        arr2d = MA.asarray(arr2d)
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"

        # check if data is balanced
        rgLens = Numeric.array(map(lambda x: len(x), replicaGroupInd), Numeric.Int)
        isBalanced = Numeric.add.reduce((1.*rgLens/rgLens[0]) == Numeric.ones((len(rgLens)))) / len(rgLens) == 1

        # check for empty cells, raise exception if empty cells exist and addInteraction=1
        if addInteraction:
            ax1Ind = Numeric.concatenate(([0], Numeric.add.accumulate(rgLens)))
            for idx in range(rgLens.shape[0]):
                if Numeric.add.reduce(MA.count(arr2d[:,ax1Ind[idx]:ax1Ind[idx+1]],1) == 0) > 0:
                    raise ValueError, "arr2d has empty cells, cannot evaluate the interaction between factors, try addInteraction=0"

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d) # [a0b0, a0b1, .., a1b0, ...]
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = numpyExtn.compressIndices(self.y)

        # dummy variables
        self.dummyA = self.getDummyEE(arr2d.shape[0], 1)
        self.dummyB = self.getDummyEE(len(replicaGroupInd), rgLens)
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

        self.FAprob = scipy.stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        self.FBprob = scipy.stats.fprob(self.dummyB.shape[1], LR_treat.DFres, self.FB)
        if addInteraction:
            self.FABprob = scipy.stats.fprob(self.dummyAB.shape[1], LR_treat.DFres, self.FAB)

        # store variables needed for posthoc tests
        self._arr2d = arr2d
        self._replicaGroupInd = replicaGroupInd
        self._LR_treat_MSres = LR_treat.MSres
        self._LR_treat_DFres = LR_treat.DFres

       
    def posthoc_tstat_B(self):
        """Fisher LSD: Compares all levels of factor B and returns (b,b) matrix of t-stat values and a matrix of corresponding p-values.
        Use when there is no significant interaction.
        Warning: Statistica computes it WRONG!
        """
        b = len(self._replicaGroupInd)
        # precompute averages
        x_avrg = Numeric.zeros((b,), Numeric.Float)         # mean of cell means over factor A
        sum_count_x = Numeric.zeros((b,), Numeric.Float)
        for idx, replicaInd in enumerate(self._replicaGroupInd):
            xi = self._arr2d.take(replicaInd, 1)
            x_avrg[idx] = MA.average(MA.average(xi,1))                  # first average over replicas to obtain cell mean, then over factor A
            sum_count_x[idx] = Numeric.add.reduce(1./MA.count(xi,1))    # first count the number of measurements within cells, then sum inverses over factor A
            ##x_avrg[idx] = MA.average(MA.ravel(xi))                    # this is how Statistica computes it, but it is WRONG!
            ##sum_count_x[idx] = 1./MA.count(xi)                        # this is how Statistica computes it, but it is WRONG!
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (x_avrg.shape[0], x_avrg.shape[0]))
        sum_count_x_2 = Numeric.resize(sum_count_x, (sum_count_x.shape[0],sum_count_x.shape[0]))
        tstat = (MA.transpose(x_avrg_2) - x_avrg_2) * self._arr2d.shape[0] / Numeric.sqrt(self._LR_treat_MSres * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))
        ##tstat = (MA.transpose(x_avrg_2) - x_avrg_2) / Numeric.sqrt(self._LR_treat_MSres * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))     # this is how Statistica computes it, but it is WRONG!
        tprob = numpyExtn.triangularPut(scipy.stats.abetai(0.5*self._LR_treat_DFres, 0.5, float(self._LR_treat_DFres) / (self._LR_treat_DFres + numpyExtn.triangularGet(tstat**2))),1,1)
        return tstat, tprob


    def isSignif_holm_B(self, alpha):
        """Conduct Holm step-down sequential multiple comparison on factor B.
        Returns (b,b) matrix indicating significant differences between levels of factor B and the direction of changes [-1|0|1].
        """
        tstat, pvals = self.posthoc_tstat_B()
        pvals = numpyExtn.triangularGet(pvals)
        sortInd = Numeric.argsort(pvals)
        k = pvals.shape[0]
        isSignif = -1*Numeric.ones((k,), Numeric.Int)
        for j in range(k):
            isSignif[sortInd[j]] = pvals[sortInd[j]] < float(alpha) / (k-j)
        return numpyExtn.triangularPut(isSignif, upper=1, lower=1) * (MA.greater(tstat,0) - MA.less(tstat,0))


    def posthoc_tstat_BbyA(self):
        """Fisher LSD: Compares all levels of factor B by each level of factor A separately.
        Returns: (b,b,a) matrix of t-stat values
                 (b,b,a) matrix of corresponding p-values
                 (b,b) matrix of geometric means of p-values over all levels of factor A.
        Use when there is a significant interaction and factor A is of no interest.
        """
        a = self._arr2d.shape[0]
        b = len(self._replicaGroupInd)
        # precompute averages
        x_avrg = Numeric.zeros((a,b), Numeric.Float)
        sum_count_x = Numeric.zeros((a,b), Numeric.Float)
        for idx, replicaInd in enumerate(self._replicaGroupInd):
            xi = self._arr2d.take(replicaInd, 1)
            x_avrg[:,idx] = MA.average(xi, 1)
            sum_count_x[:,idx] = 1. / MA.count(xi, 1)
        
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (b,a,b))
        sum_count_x_2 = Numeric.resize(sum_count_x, (b,a,b))

        tstat = (MA.transpose(x_avrg_2, (2,1,0)) - x_avrg_2) / Numeric.sqrt(self._LR_treat_MSres * (Numeric.transpose(sum_count_x_2, (2,1,0)) + sum_count_x_2))
        # get p-values for each level of factor a separately (axis 1)
        tprob = MA.array(-1*Numeric.ones(tstat.shape, Numeric.Float), mask=Numeric.ones(tstat.shape))
        for idx1 in range(tprob.shape[1]):
            tprob[:,idx1,:] = numpyExtn.triangularPut(scipy.stats.abetai(0.5*self._LR_treat_DFres, 0.5, float(self._LR_treat_DFres) / (self._LR_treat_DFres + numpyExtn.triangularGet(tstat[:,idx1,:]**2))),1,1)
        return tstat, tprob


    def posthoc_anova_B_AB(self):
        """Conduct ANOVA for each pair of levels of factor B to detect interaction effect.
        Returns for each pair of factor B levels:
            (b,b) matrix of F-stat values for effect of factor B 
            (b,b) matrix of p-values for effect of factor B
            (b,b) matrix of F-stat values for interaction effect
            (b,b) matrix of p-values for interaction effect
        """
        b = len(self._replicaGroupInd)
        FB = MA.masked * MA.ones((b,b), Numeric.Float)
        FBprob = MA.masked * MA.ones((b,b), Numeric.Float)
        FAB = MA.masked * MA.ones((b,b), Numeric.Float)
        FABprob = MA.masked * MA.ones((b,b), Numeric.Float)
        for i in range(b):
            for j in range(i+1,b):
                rglSub = [self._replicaGroupInd[i], self._replicaGroupInd[j]]
                takeInd = self._replicaGroupInd[i] + self._replicaGroupInd[j]
                an = Anova2wayLR(self._arr2d.take(takeInd, 1), rglSub, addInteraction=1)
                FB[i,j] = FB[j,i] = an.FB
                FBprob[i,j] = FBprob[j,i] = an.FBprob
                FAB[i,j] = FAB[j,i] = an.FAB
                FABprob[i,j] = FABprob[j,i] = an.FABprob
        return FB, FBprob, FAB, FABprob


class Anova2wayLR(Anova2wayLRBase):
    """2 way ANOVA with multiple observations per cell for unbalanced designs using multivariate linear regression.
    Multiple observations given at axis 1 together with the list of group lenghts, e.g. [3,2].
    Supports balanced and unbalanced designs, i.e. different number of observations per cell and missing (masked) data.
    TODO: account for empty cells (this implementation results in singular design matrix)
    """

    def __init__(self, arr2d, groupLens, addInteraction=0, allowReductA=True, allowReductB=False):
        """arr2d: 2D masked array where replicas are given at axis 1 together with lenghts of the groups (multiple observations per cell);
        addInteraction: include / exclude the interaction between factors;
        allowReduct[A|B]: whether to allow or not the reduction of levels of factor A (B).
        """
        arr2d = MA.asarray(arr2d)
        groupLens = Numeric.array(groupLens)    # make a copy for the case if subjects or levels of factor B are removed
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"

        # check if there exist factor-levels with all values missing (for A and B)
        # if so, reduce the corresponding DF and fix the number of dummy variables
        missIndA = Numeric.compress(MA.count(arr2d, 1) == 0, Numeric.arange(arr2d.shape[0]))
        if missIndA.shape[0] > 0:
            if allowReductA:
                print "Warning: removig factor A level(s) %s" % str(missIndA.tolist())
                takeIndA = Numeric.compress(MA.count(arr2d, 1) != 0, Numeric.arange(arr2d.shape[0]))
                arr2d = arr2d.take(takeIndA, 0)
            else:
                self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: the following level(s) of factor A have no values: %s" % str(missIndA.tolist()))
                return
        missIndSubj = Numeric.compress(MA.count(arr2d, 0) == 0, Numeric.arange(arr2d.shape[1]))
        if missIndSubj.shape[0] > 0:
            takeIndSubj = Numeric.compress(MA.count(arr2d, 0) != 0, Numeric.arange(arr2d.shape[1]))
            # fix groupLens
##            mapSubj2Group = -1*Numeric.ones(arr2d.shape[1])
##            groupLensAcc0 = [0] + Numeric.add.accumulate(groupLens).tolist()
##            for i in range(len(groupLensAcc0) - 1):
##                mapSubj2Group[groupLensAcc0[i]:groupLensAcc0[i+1]] = i
            mapSubj2Group = Numeric.repeat(range(len(groupLens)), groupLens)
            for subjIdx in missIndSubj:
                groupLens[mapSubj2Group[subjIdx]] -= 1
            # remove data columns that are missing all the values
            arr2d = arr2d.take(takeIndSubj, 1)
        # fix number of factor B levels (if the length of any group became 0 due to removed subjects)
        missIndB = Numeric.compress(groupLens <= 0, Numeric.arange(groupLens.shape[0]))
        if missIndB.shape[0] > 0:
            if allowReductB:
                print "Warning: removig factor B level(s) %s" % str(missIndB)
                takeIndB = Numeric.compress(groupLens > 0, Numeric.arange(groupLens.shape[0]))
                groupLens = Numeric.take(groupLens, takeIndB)
                # arr2d has already been taken care of by removing the subjects without observations
            else:
                self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: the following level(s) of factor B have no values: %s" % str(missIndB))
                return
        # remove levels of factor B with single observation per cell
        grpIndSingle = Numeric.compress(Numeric.equal(groupLens,1), Numeric.arange(groupLens.shape[0]))
        if grpIndSingle.shape[0] > 0:
            if allowReductB:
                # fix arr2d
                print "Warning: removing factor B level(s) with single observation: %s" % str(grpIndSingle)
                groupLensAcc = [0] + list(Numeric.add.accumulate(groupLens))
                arr2dRemInd = map(lambda i: groupLensAcc[i], grpIndSingle)
                arr2dTakeInd = Numeric.ones(arr2d.shape[1])
                Numeric.put(arr2dTakeInd, arr2dRemInd, 0)
                arr2d = MA.compress(arr2dTakeInd, arr2d, 1)
                # fix groupLens
                groupLens = Numeric.compress(Numeric.not_equal(groupLens,1), groupLens)
            else:
                self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: the following level(s) of factor B have single observations: %s" % str(grpIndSingle))
                return                
        # check that there exist at least 2 levels of factors A and B
        # and that there is at least single duplicated observation
        if arr2d.shape[0] < 2 or len(groupLens) < 2: ##or MA.count(arr2d) <= arr2d.shape[0] * len(groupLens):
            self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: not enough levels of factor A (%i) or factor B (%i)" % (arr2d.shape[0],len(groupLens)))
            return

        # check if data is balanced
        isBalanced = (1. * Numeric.add.reduce((1.*groupLens/groupLens[0]) == Numeric.ones(groupLens.shape[0])) / groupLens.shape[0]) == 1

        # check for empty cells; if exist and addInteraction=1, switch to no interaction and continue
        if addInteraction:
            ax1Ind = Numeric.concatenate(([0], Numeric.add.accumulate(groupLens)))
            for idx in range(groupLens.shape[0]):
                if Numeric.add.reduce(MA.count(arr2d[:,ax1Ind[idx]:ax1Ind[idx+1]],1) == 0) > 0:
##                    raise ValueError, "arr2d has empty cells, cannot evaluate the interaction between factors, try addInteraction=0"
                    print ValueError, "Warning: data has empty cells, cannot evaluate the interaction between factors, continuing without interaction"
                    addInteraction=0
                    # store member varibales that are stored only if addInteraction==1
                    self.FAB = 0
                    self.FABprob = 1
                    self.addInteractionFailed = True

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d) # [a0b0, a0b1, .., a1b0, ...]
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = numpyExtn.compressIndices(self.y)

        # dummy variables
        self.dummyA = self.getDummyEE(arr2d.shape[0], 1)
        self.dummyB = self.getDummyEE(len(groupLens), groupLens)
        self.dummyA, self.dummyB = self.getDummiesJoin(self.dummyA, self.dummyB)
        if addInteraction:
            self.dummyAB = self.getDummyInteraction(self.dummyA, self.dummyB)
        else:
            self.dummyAB = Numeric.zeros((self.dummyA.shape[0],0))
        self.dummyA = Numeric.take(self.dummyA, takeInd, 0)
        self.dummyB = Numeric.take(self.dummyB, takeInd, 0)
        self.dummyAB = Numeric.take(self.dummyAB, takeInd, 0)

        # check that there is enough observations to allow for variability (DFres > 0)
        if self.dummyA.shape[0] - sum([self.dummyA.shape[1], self.dummyB.shape[1], self.dummyAB.shape[1]]) - 1 <= 0:
            self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: not values (dummy condition)")
            return

        # run analysis
        try:
##            if addInteraction:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB, self.dummyAB), 1), self.y)                
##            else:
##                LR_treat = MultLinReg(Numeric.concatenate((self.dummyA, self.dummyB),1), self.y)
                
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
        except:
            self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: miltivariate linear regression failed due to lack of variability")
            return

        try:
            self.FAprob = scipy.stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        except ValueError:
            self.FAprob = 1
        try:
            self.FBprob = scipy.stats.fprob(self.dummyB.shape[1], LR_treat.DFres, self.FB)
        except ValueError:
            self.FBprob = 1
        self.ps = [self.FAprob, self.FBprob]
        if addInteraction:
            try:
                self.FABprob = scipy.stats.fprob(self.dummyAB.shape[1], LR_treat.DFres, self.FAB)
            except ValueError:
                self.FABprob = 1
            self.ps.append(self.FABprob)

        # store variables needed for posthoc tests
        self._MSpool = LR_treat.MSres
        self._DFpool = LR_treat.MSreg
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens


    def _giveUp(self, arr2d, groupLens, addInteraction, output=None):
        if output != None:
            print output
        self.FA = 0
        self.FAprob = 1
        self.FB = 0
        self.FBprob = 1
        self.ps = [self.FAprob, self.FBprob]
        if addInteraction:
            self.FAB = 0
            self.FABprob = 1
            self.ps.append(self.FABprob)
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens
        


class AnovaRM11LR(AnovaLRBase):
    """1 way REPEATED MEASURES ANOVA using multivariate linear regression.
    Axis 0 correspond to subjects, axis 1 correspond to different factor-levels.
    Does NOT support missing (masked) data.
    """
    
    def __init__(self, arr2d):
        """arr2d: 2D masked array: [[x1,x2], [y1,y2], [z1,z2]]
        where xi,yi,zi correspond to measurements of three different subjects (x,y,z) under the i-th factor-level.
        """
        arr2d = MA.asarray(arr2d)
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d)
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = numpyExtn.compressIndices(self.y)

        # dummy variables for coding subjects and factor-levels
        self.dummySubj = self.getDummyEE(arr2d.shape[0], 1)
        self.dummyA = self.getDummyEE(arr2d.shape[1], 1)
        self.dummySubj, self.dummyA = self.getDummiesJoin(self.dummySubj, self.dummyA)
        self.dummySubj = Numeric.take(self.dummySubj, takeInd, 0)
        self.dummyA = Numeric.take(self.dummyA, takeInd, 0)

        # run analysis (better: only 2 LR instead of 3)
        try:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummySubj, self.dummyA), 1), self.y)
            if noMissing:
                LR_A = MultLinReg(self.dummyA, self.y)
                self.FA = LR_A.MSreg / LR_treat.MSres
            else:
                print "WARNING: missing data with 1-way repeated measures ANOVA"
                LR_S = MultLinReg(self.dummySubj, self.y)
                self.FA = (LR_treat.SSreg - LR_S.SSreg) / self.dummyA.shape[1] / LR_treat.MSres
            self.FAprob = scipy.stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        except ZeroDivisionError:
            self.FA = 0
            self.FAprob = 1


class _AnovaRM12LR_bug_missing_factor_leves(AnovaLRBase):
    """2 way ANOVA with REPEATED MEASURES on factor A using multivariate linear regression.
    Axis 0 correspond to different levels of factor A, subjects given at axis 1, subjects nested inside factor B according to the given group indices.
    Factor A is a within-subject effect, factor B is a between-ubject effect.
    Example:
        arr2d = [[a0b0 a0b0 a0b1 a0b1 a0b1], [a1b0 a1b0 a1b1 a1b1 a1b1]]: 2 levels of factor A, 2 levels of factor B, 5 subjects
        groupInd = [[0,1],[2,3,4]]
    Supports balanced and unbalanced designs, i.e. different number of observations per cell and missing (masked) data.
    TODO: fix the F and DF for factor B when there are missing values (see Glantz, pp. 462-5)
    """

    def __init__(self, arr2d, groupInd, addInteraction=0):
        """arr2d: 2D masked array: [[a0b0 a0b0 a0b1 a0b1 a0b1], [a1b0 a1b0 a1b1 a1b1 a1b1]]
        groupInd: axis 1 indices that belong to the same level of factor B: [[0,1],[2,3,4]]
        addInteraction: include / exclude the interaction between factors A and B
        """
        arr2d = MA.asarray(arr2d)
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"

        # remove missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d) # [a0b0, a0b1, .., a1b0, ...]
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = numpyExtn.compressIndices(self.y)

        # check if data is balanced
        rgLens = Numeric.array(map(lambda x: len(x), groupInd), Numeric.Int)
        assert Numeric.add.reduce(rgLens) == arr2d.shape[1], "arr2d.shape[1] != sum_of_len(groupInd)"
        isBalanced = Numeric.add.reduce((1.*rgLens/rgLens[0]) == Numeric.ones((len(rgLens)))) / len(rgLens) == 1

        # dummy variables
        self.dummyA = self.getDummyEE(arr2d.shape[0], 1)
        self.dummyB = self.getDummyEE(len(groupInd), rgLens)
        self.dummySubjInB = self.getSubjectDummyEE(rgLens, 1)
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

        self.FBprob = scipy.stats.fprob(self.dummyB.shape[1], self.dummySubjInB.shape[1], self.FB)
        self.FAprob = scipy.stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        if addInteraction:
            self.FABprob = scipy.stats.fprob(self.dummyAB.shape[1], LR_treat.DFres, self.FAB)

        # store variables needed for posthoc tests
        self._arr2d = arr2d
        self._groupInd = groupInd
        self._LR_treat_MSres = LR_treat.MSres
        self._LR_treat_DFres = LR_treat.DFres
        # estimate the variance by pooling residual variance for factors A and B
        MS_res = LR_treat.MSres
        a = self.dummyA.shape[1] + 1    # number of levels of factor A
        DF_SubjInB = self.dummySubjInB.shape[1] # = LR_SubjInB.DFreg
        DF_res = LR_treat.DFres
        self._MSpool = (MS_SubjInB + (a-1)*MS_res) / a
        self._DFpool = ((MS_SubjInB + (a-1)*MS_res)**2) / ((MS_SubjInB**2 / DF_SubjInB) + (((a-1)*MS_res)**2 / DF_res))


    def posthoc_tstat_B(self):
        """Fisher LSD: Compares all levels of factor B and returns (b,b) matrix of t-stat values and a matrix of corresponding p-values.
        Use when there is no significant interaction.
        Warning: Statistica computes it WRONG!
        """
        b = len(self._groupInd)
        # precompute averages
        x_avrg = Numeric.zeros((b,), Numeric.Float)
        sum_count_x = Numeric.zeros((b,), Numeric.Float)
        for idx, replicaInd in enumerate(self._groupInd):
            xi = self._arr2d.take(replicaInd, 1)
            x_avrg[idx] = MA.average(MA.average(xi,1))                  # first average over replicas to obtain cell mean, then over factor A
            sum_count_x[idx] = Numeric.add.reduce(1./MA.count(xi,1))    # first count the number of measurements within cells, then sum inverses over factor A
            ##x_avrg[idx] = MA.average(MA.ravel(xi))                    # this is how Statistica computes it, but it is WRONG!
            ##sum_count_x[idx] = 1./MA.count(xi)                        # this is how Statistica computes it, but it is WRONG!
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (x_avrg.shape[0], x_avrg.shape[0]))
        sum_count_x_2 = Numeric.resize(sum_count_x, (sum_count_x.shape[0],sum_count_x.shape[0]))
        tstat = (MA.transpose(x_avrg_2) - x_avrg_2) * self._arr2d.shape[0] / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))
        ##tstat = (MA.transpose(x_avrg_2) - x_avrg_2) / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2) + sum_count_x_2))     # this is how Statistica computes it, but it is WRONG!
        tprob = numpyExtn.triangularPut(scipy.stats.abetai(0.5*self._DFpool, 0.5, float(self._DFpool) / (self._DFpool + numpyExtn.triangularGet(tstat**2))),1,1)
        return tstat, tprob
        
            
    def posthoc_tstat_BbyA(self):
        """Fisher LSD: Compares all levels of factor B by each level of factor A separately.
        Returns: (b,b,a) matrix of t-stat values
                 (b,b,a) matrix of corresponding p-values
                 (b,b) matrix of geometric means of p-values over all levels of factor A.
        Use when there is a significant interaction and factor A is of no interest.
        """
        a = self._arr2d.shape[0]
        b = len(self._groupInd)
        # precompute averages
        x_avrg = Numeric.zeros((a,b), Numeric.Float)
        sum_count_x = Numeric.zeros((a,b), Numeric.Float)
        for idx, replicaInd in enumerate(self._groupInd):
            xi = self._arr2d.take(replicaInd, 1)
            x_avrg[:,idx] = MA.average(xi, 1)
            sum_count_x[:,idx] = 1. / MA.count(xi, 1)
        
        # t-statistics
        x_avrg_2 = MA.resize(x_avrg, (b,a,b))
        sum_count_x_2 = Numeric.resize(sum_count_x, (b,a,b))

        tstat = (MA.transpose(x_avrg_2, (2,1,0)) - x_avrg_2) / Numeric.sqrt(self._MSpool * (Numeric.transpose(sum_count_x_2, (2,1,0)) + sum_count_x_2))
        # get p-values for each level of factor a separately (axis 1)
        tprob = MA.array(-1*Numeric.ones(tstat.shape, Numeric.Float), mask=Numeric.ones(tstat.shape))
        for idx1 in range(tprob.shape[1]):
            tprob[:,idx1,:] = numpyExtn.triangularPut(scipy.stats.abetai(0.5*self._DFpool, 0.5, float(self._DFpool) / (self._DFpool + numpyExtn.triangularGet(tstat[:,idx1,:]**2))),1,1)
        return tstat, tprob

    def posthoc_anova_B_AB(self):
        """Conduct ANOVA for each pair of levels of factor B to detect interaction effect.
        Returns for each pair of factor B levels:
            (b,b) matrix of F-stat values for effect of factor B 
            (b,b) matrix of p-values for effect of factor B
            (b,b) matrix of F-stat values for interaction effect
            (b,b) matrix of p-values for interaction effect
        """
        b = len(self._groupInd)
        FB = MA.masked * MA.ones((b,b), Numeric.Float)
        FBprob = MA.masked * MA.ones((b,b), Numeric.Float)
        FAB = MA.masked * MA.ones((b,b), Numeric.Float)
        FABprob = MA.masked * MA.ones((b,b), Numeric.Float)
        for i in range(b):
            for j in range(i+1,b):
                rglSub = [self._groupInd[i], self._groupInd[j]]
                takeInd = self._groupInd[i] + self._groupInd[j]
                an = AnovaRM12LR(self._arr2d.take(takeInd, 1), rglSub, addInteraction=1)
                FB[i,j] = FB[j,i] = an.FB
                FBprob[i,j] = FBprob[j,i] = an.FBprob
                FAB[i,j] = FAB[j,i] = an.FAB
                FABprob[i,j] = FABprob[j,i] = an.FABprob
        return FB, FBprob, FAB, FABprob
    

class AnovaRM12LR(Anova2wayLRBase):
    """2 way ANOVA with REPEATED MEASURES on factor A using multivariate linear regression.
    Axis 0 correspond to different levels of factor A, subjects given at axis 1, subjects nested inside factor B according to the given group lenghts.
    Factor A is a within-subject effect, factor B is a between-subject effect.
    Example:
        arr2d = [[a0b0 a0b0 a0b1 a0b1 a0b1], [a1b0 a1b0 a1b1 a1b1 a1b1]]: 2 levels of factor A, 2 levels of factor B, 5 subjects
        groupLens = [2,3]
    Supports balanced and unbalanced designs, i.e. different number of observations per cell and missing (masked) data.
    TODO: fix the F and DF for factor B when there are missing values (see Glantz, pp. 462-5)
    """

    def __init__(self, arr2d, groupLens, addInteraction=0, allowReductA=True, allowReductB=False):
        """arr2d: 2D masked array: [[a0b0 a0b0 a0b1 a0b1 a0b1], [a1b0 a1b0 a1b1 a1b1 a1b1]];
        groupLens: lenghts of axis 1 indices that belong to the same level of factor B, e.g. [2,3];
        addInteraction: include / exclude the interaction between factors A and B;
        allowReduct[A|B]: whether to allow or not the reduction of levels of factor A (B).
        """
        arr2d = MA.asarray(arr2d)
        groupLens = Numeric.array(groupLens)    # make a copy for the case if subjects or levels of factor B are removed
        assert len(arr2d.shape) == 2, "len(arr2d.shape) != 2"
        assert arr2d.shape[1] == Numeric.add.reduce(groupLens), "arr2d.shape[1] != Numeric.add.reduce(groupLens)"

        # check if there exist factor-levels with all values missing (for A, B and Subj)
        # if so, reduce the corresponding DF and fix the number of dummy variables
        missIndA = Numeric.compress(MA.count(arr2d, 1) == 0, Numeric.arange(arr2d.shape[0]))
        if missIndA.shape[0] > 0:
            if allowReductA:
                print "Warning: removig factor A level(s) %s" % str(missIndA.tolist())
                takeIndA = Numeric.compress(MA.count(arr2d, 1) != 0, Numeric.arange(arr2d.shape[0]))
                arr2d = arr2d.take(takeIndA, 0)
            else:
                self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: the following level(s) of factor A have no values: %s" % str(missIndA.tolist()))
                return
        missIndSubj = Numeric.compress(MA.count(arr2d, 0) == 0, Numeric.arange(arr2d.shape[1]))
        if missIndSubj.shape[0] > 0:
            print "Warning: removig subject(s) %s" % str(missIndSubj.tolist())
            takeIndSubj = Numeric.compress(MA.count(arr2d, 0) != 0, Numeric.arange(arr2d.shape[1]))
            # fix groupLens
##            mapSubj2Group = -1*Numeric.ones(arr2d.shape[1])
##            groupLensAcc0 = [0] + Numeric.add.accumulate(groupLens).tolist()
##            for i in range(len(groupLensAcc0) - 1):
##                mapSubj2Group[groupLensAcc0[i]:groupLensAcc0[i+1]] = i
            mapSubj2Group = Numeric.repeat(range(len(groupLens)), groupLens)
            for subjIdx in missIndSubj:
                groupLens[mapSubj2Group[subjIdx]] -= 1
            # remove data columns that are missing all the values
            arr2d = arr2d.take(takeIndSubj, 1)
        # fix number of factor B levels (if the length of any group became 0 due to removed subjects)
        missIndB = Numeric.compress(groupLens <= 0, Numeric.arange(groupLens.shape[0]))
        if missIndB.shape[0] > 0:
            if allowReductB:
                print "Warning: removig factor B level(s) %s" % str(missIndB.tolist())
                takeIndB = Numeric.compress(groupLens > 0, Numeric.arange(groupLens.shape[0]))
                groupLens = Numeric.take(groupLens, takeIndB)
                # arr2d has already been taken care of by removing the subjects without observations
            else:
                self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: the following level(s) of factor B have no values: %s" % str(missIndB))
                return
        # remove levels of factor B with single observation per cell
        grpIndSingle = Numeric.compress(Numeric.equal(groupLens,1), Numeric.arange(groupLens.shape[0]))
        if grpIndSingle.shape[0] > 0:
            if allowReductB:
                # fix arr2d
                print "Warning: removing factor B level(s) with single observation: %s" % str(grpIndSingle)
                groupLensAcc = [0] + list(Numeric.add.accumulate(groupLens))
                arr2dRemInd = map(lambda i: groupLensAcc[i], grpIndSingle)
                arr2dTakeInd = Numeric.ones(arr2d.shape[1])
                Numeric.put(arr2dTakeInd, arr2dRemInd, 0)
                arr2d = MA.compress(arr2dTakeInd, arr2d, 1)
                # fix groupLens
                groupLens = Numeric.compress(Numeric.not_equal(groupLens,1), groupLens)
            else:
                self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: the following level(s) of factor B have single observations: %s" % str(grpIndSingle))
                return                
        # check that there exist at least 2 levels of factors A and B
        # and that there is at least single duplicated observation
        if arr2d.shape[0] < 2 or len(groupLens) < 2: ## or MA.count(arr2d) <= arr2d.shape[0] * len(groupLens):
            self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: not enough levels of factor A (%i) or factor B (%i)" % (arr2d.shape[0],len(groupLens)))
            return

        # remove other missing observations and the corresponding dummy variable codes
        self.y = MA.ravel(arr2d) # [a0b0, a0b1, .., a1b0, ...]
        noMissing = MA.count(self.y) == self.y.shape[0]
        self.y, takeInd = numpyExtn.compressIndices(self.y)

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
        else:
            self.dummyAB = Numeric.zeros((self.dummyA.shape[0],0))
        self.dummyA = Numeric.take(self.dummyA, takeInd, 0)
        self.dummyB = Numeric.take(self.dummyB, takeInd, 0)
        self.dummySubjInB = Numeric.take(self.dummySubjInB, takeInd, 0)

        # check that there is enough observations to allow for variability (DFres > 0)
        if self.dummyA.shape[0] - sum([self.dummySubjInB.shape[1], self.dummyA.shape[1], self.dummyB.shape[1], self.dummyAB.shape[1]]) - 1 <= 0:
            self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: not values (dummy condition)")
            return

        # run analysis
        try:
##            if addInteraction:
            LR_treat = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA, self.dummyB, self.dummyAB),1), self.y)
##            else:
##                LR_treat = MultLinReg(Numeric.concatenate((self.dummySubjInB, self.dummyA, self.dummyB),1), self.y)            

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
        except:
            self._giveUp(arr2d, groupLens, addInteraction, output="Giving up ANOVA: miltivariate linear regression failed, probably due to lack of variability")
            return
        
        try:
            self.FBprob = scipy.stats.fprob(self.dummyB.shape[1], self.dummySubjInB.shape[1], self.FB)
        except ValueError:
            self.FBprob = 1
        try:
            self.FAprob = scipy.stats.fprob(self.dummyA.shape[1], LR_treat.DFres, self.FA)
        except ValueError:
            self.FAprob = 1
        self.ps = [self.FAprob, self.FBprob]
        if addInteraction:
            try:
                self.FABprob = scipy.stats.fprob(self.dummyAB.shape[1], LR_treat.DFres, self.FAB)
            except ValueError:
                self.FABprob = 1
            self.ps.append(self.FABprob)

        # estimate the variance and degrees of freedom by pooling residual variance for factors A and B
        MS_res = LR_treat.MSres
        a = self.dummyA.shape[1] + 1    # number of levels of factor A
        DF_SubjInB = self.dummySubjInB.shape[1] # = LR_SubjInB.DFreg
        DF_res = LR_treat.DFres
        self._MSpool = (MS_SubjInB + (a-1)*MS_res) / a
        self._DFpool = ((MS_SubjInB + (a-1)*MS_res)**2) / ((MS_SubjInB**2 / DF_SubjInB) + (((a-1)*MS_res)**2 / DF_res))

        # store some variables
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens


    def _giveUp(self, arr2d, groupLens, addInteraction, output=None):
        if output != None:
            print output
        self.FA = 0
        self.FAprob = 1
        self.FB = 0
        self.FBprob = 1
        self.ps = [self.FAprob, self.FBprob]
        if addInteraction:
            self.FAB = 0
            self.FABprob = 1
            self.ps.append(self.FABprob)
        self._addInteraction = addInteraction
        self._arr2d = arr2d
        self._groupLens = groupLens


#######################################################################################
## Traditional ANOVA
#######################################################################################

class Anova2way:
    """2 way ANOVA with multiple observations per cell for balanced designs.
    Multiple observations given at axis 1 together with the number of repeats.
    """

    def __init__(self, arr2d, numRepeats):
        """Input: 2d array with multiple observations per cell ginven at axis 1.
        """
        arr2d = MA.asarray(arr2d)
        assert math.modf(arr2d.shape[1] / numRepeats)[0] == 0.0, "arr2d.shape[1] / numRepeats not a whole number"
        self.a = arr2d.shape[0]
        self.b = arr2d.shape[1] / numRepeats
        self.n = numRepeats     # number of observations per cell
        self.Y = arr2d
        self._init_sums()
        self._init_SS()
        self._init_F()


    def _init_sums(self):
        """Initializes sums of Y.
        """
        self.dfA = self.a - 1
        self.dfB = self.b - 1
        self.dfAB = self.dfA * self.dfB
        self.dfTotal = self.a * self.b * self.n  - 1
        self.dfError = self.dfTotal - self.dfA - self.dfB - self.dfAB
        # YijK
        self.YijK = Numeric.zeros((self.a,self.b), Numeric.Float)
        for idx, bIdx in enumerate(range(0, self.b*self.n, self.n)):
            self.YijK[:,idx] = Numeric.add.reduce(self.Y[:,bIdx:bIdx+self.n], 1)
        # YiJK
        self.YiJK = Numeric.add.reduce(self.YijK, 1)
        # YIjK
        self.YIjK = Numeric.add.reduce(self.YijK, 0)
        # YIJK 
        self.YIJK = Numeric.add.reduce(self.YIjK, 0)

        
    def _init_SS(self):
        y2IJKdivnab = Numeric.power(self.YIJK,2) / (self.n * self.a * self.b)
        self.TSS = Numeric.add.reduce(Numeric.add.reduce(Numeric.power(self.Y,2))) - y2IJKdivnab
        self.SSA = Numeric.add.reduce(Numeric.power(self.YiJK,2)) / (self.n * self.b) - y2IJKdivnab
        self.SSB = Numeric.add.reduce(Numeric.power(self.YIjK,2)) / (self.n * self.a) - y2IJKdivnab
        self.SSAB = Numeric.add.reduce(Numeric.add.reduce(Numeric.power(self.YijK,2))) / (self.n) - y2IJKdivnab - self.SSA - self.SSB
        self.SSE = self.TSS - self.SSA - self.SSB - self.SSAB
        

    def _init_F(self):
        self.FA = (self.SSA / self.dfA) / (self.SSE / self.dfError)
        self.FAprob = scipy.stats.fprob(self.dfA, self.dfError, self.FA)
        self.FB = (self.SSB / self.dfB) / (self.SSE / self.dfError)
        self.FBprob = scipy.stats.fprob(self.dfB, self.dfError, self.FB)
        self.FAB = (self.SSAB / self.dfAB) / (self.SSE / self.dfError)
        self.FABprob = scipy.stats.fprob(self.dfAB, self.dfError, self.FAB)
        
        
class Anova3way:
    """3 way ANOVA, optional AB interaction, single observation per cell, balanced design.
    Removes factor-levels with mising observations.
    """

    def __init__(self, arr3d, addInterAB):
        """Input: 3d masked array where axis correspond factors.
        Leave out the factor-levels with masked values.
        """
        arr3d = MA.asarray(arr3d)
        assert len(arr3d.shape) == 3, "len(arr3d.shape) != 3"
        a,b,c = arr3d.shape
        self.Y = arr3d

        self.YijK = Numeric.add.reduce(self.Y, 2)
        self.YiJK = Numeric.add.reduce(self.YijK, 1)
        self.YIjK = Numeric.add.reduce(self.YijK, 0)
        self.YIJk = Numeric.add.reduce(Numeric.add.reduce(self.Y))
        self.YIJK2_abc = Numeric.power(Numeric.add.reduce(self.YIJk),2) / (a*b*c)

        self.TSS = Numeric.add.reduce(Numeric.add.reduce(Numeric.add.reduce(Numeric.power(self.Y,2)))) - self.YIJK2_abc
        self.SSA = Numeric.add.reduce(Numeric.power(self.YiJK,2)) / (b*c) - self.YIJK2_abc
        self.SSB = Numeric.add.reduce(Numeric.power(self.YIjK,2)) / (a*c) - self.YIJK2_abc
        self.SSC = Numeric.add.reduce(Numeric.power(self.YIJk,2)) / (a*b) - self.YIJK2_abc

        self.dfA = a - 1
        self.dfB = b - 1
        self.dfC = c - 1
        self.dfTotal = a * b * c  - 1

        if addInterAB:

            self.SSAB = Numeric.add.reduce(Numeric.add.reduce(Numeric.power(self.YijK,2))) / c - self.YIJK2_abc - self.SSA - self.SSB
            self.dfAB = self.dfA * self.dfB
            
            self.SSE = self.TSS - self.SSA - self.SSB - self.SSC - self.SSAB
            self.dfError = self.dfTotal - self.dfA - self.dfB - self.dfC - self.dfAB

            self.FAB = (self.SSAB / self.dfAB) / (self.SSE / self.dfError)
            self.FABprob = scipy.stats.fprob(self.dfAB, self.dfError, self.FAB)

        else:            

            self.SSE = self.TSS - self.SSA - self.SSB - self.SSC
            self.dfError = self.dfTotal - self.dfA - self.dfB - self.dfC

        self.FA = (self.SSA / self.dfA) / (self.SSE / self.dfError)
        self.FAprob = scipy.stats.fprob(self.dfA, self.dfError, self.FA)
        self.FB = (self.SSB / self.dfB) / (self.SSE / self.dfError)
        self.FBprob = scipy.stats.fprob(self.dfB, self.dfError, self.FB)
        self.FC = (self.SSC / self.dfC) / (self.SSE / self.dfError)
        self.FCprob = scipy.stats.fprob(self.dfC, self.dfError, self.FC)


#######################################################################################
## Multivariate Linear Regression
#######################################################################################


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
##        self.MSres = self.SSres / self.DFres                    # equals to square std. error: s^2_{y|x}
##        if self.MSres == 0.0:
##            print "Warning MultLinReg: MSres equals 0, replaced by 1e-20"
##            self.MSres = 1e-20
        if self.DFres > 0:
            self.MSres = self.SSres / self.DFres                    # equals to square std. error: s^2_{y|x}
        else:
            print "Warning MultLinReg: SSres=%f, DFres=%f" % (self.SSres, self.DFres)
            self.MSres = 0
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
        self.b_tprob = scipy.stats.abetai(0.5*self.DFres, 0.5, float(self.DFres)/(self.DFres + self.b_t**2))
        # correlation between coefficients (b's)
        # ERROR in Glantz et al.: self.Corr = Numeric.sqrt(self.VarCovar) / Numeric.sqrt(Numeric.dot(self.b_se[:,Numeric.NewAxis],self.b_se[Numeric.NewAxis,:]))
        self.Corr = self.VarCovar / Numeric.dot(self.b_se[:,Numeric.NewAxis],self.b_se[Numeric.NewAxis,:])
        self.F = self.MSreg / self.MSres        # overall goodness of fit 
        self.Fprob = scipy.stats.fprob(self.DFreg, self.DFres, self.F)
        self.R2 = self.SSreg/self.SStot         # coef. of determination: fraction of variance explained by regression plane
        self.R2_adj = 1-(self.MSres/self.MStot) # adjusted coef. of determination (unbiased estimator of population R2; larget the better)
        self.R = math.sqrt(self.R2)             # multiple correlation coefficeint: how tight the data points cluster around the regression plane

        

#######################################################################################
## FDR and q-value
#######################################################################################

def qVals(pVals, estimatePi0=False, verbose=False):
    """Returns q-values calculated from given p-values.
    Input:  array of p-values.
            estimate pi0: False (pi0==1)
                          "spline": cubic spline fit on pi0(lambda) curve, estimete at lambda=1
                          "loess": loess fit on pi0(lambda) curve, estimate at lambda=1
    Reference:  Storey et al. (2003) Statistical significance for genomewide studies.
                PNAS 100(16), pp. 9440-5.
    Changes: loess fit instead of cubic spline
    Added:  handles missing values -> it compresses it out
            lambda values depend on pValues:
                min(pVals) <= lmbd < max(pVals)
                len(lmbd) == len(pVals)
    TODO: remove calculation of flmbd in all points except for lmbd=1
    """
    assert estimatePi0 in [False, "spline", "loess"]
    pValsMA = MA.asarray(pVals)
    qVals = MA.ones(pValsMA.shape, Numeric.Float) * MA.masked
    putInd = Numeric.compress(Numeric.logical_not(MA.getmaskarray(pValsMA)), Numeric.arange(pValsMA.shape[0]))
    pValsMA = pValsMA.compressed()
    minp = min(pValsMA); maxp = max(pValsMA)
    m = pValsMA.shape[0]
    x = Numeric.arange(0,1.01,0.01)
    if estimatePi0 == "spline":
        import scipy.interpolate
        lmbd = Numeric.arange(0,0.96,0.01,Numeric.Float)
        pi0lmbd = map(lambda l: Numeric.add.reduce(Numeric.greater(pValsMA,l)) / (m*(1-l)), lmbd)
        splineRep = scipy.interpolate.splrep(lmbd, pi0lmbd, k=3, s=0.01)   # spline representation: (knots, coefficeints, degree)
        flmbd = scipy.interpolate.splev(x, splineRep, der=0)
        pi0 = flmbd[-1]
        if pi0 <= 0:
            print "Warning: pi0<=0, trying smaller lambda..."
            while pi0 <= 0:
                lmbd = lmbd[:-1]
                pi0lmbd = map(lambda l: Numeric.add.reduce(Numeric.greater(pValsMA,l)) / (m*(1-l)), lmbd)
                splineRep = scipy.interpolate.splrep(lmbd, pi0lmbd, k=3, s=0.1)   # spline representation: (knots, coefficeints, degree)
                flmbd = scipy.interpolate.splev(x, splineRep, der=0)
                pi0 = flmbd[-1]
    elif estimatePi0 == "loess":
        lmbd = Numeric.arange(0,1.0,0.01,Numeric.Float)
        pi0lmbd = map(lambda l: Numeric.add.reduce(Numeric.greater(pValsMA,l)) / (m*(1-l)), lmbd)
        flmbd = Numeric.asarray(statc.loess(zip(lmbd, pi0lmbd), list(x), 0.4))[:,1]
        pi0 = flmbd[-1]
        if pi0 <= 0:
            print "Warning: pi0<=0, trying smaller lambda..."
            while pi0 <= 0:
                lmbd = lmbd[:-1]
                pi0lmbd = map(lambda l: Numeric.add.reduce(Numeric.greater(pValsMA,l)) / (m*(1-l)), lmbd)
                flmbd = Numeric.asarray(statc.loess(zip(lmbd, pi0lmbd), list(x), 0.4))[:,1]
                pi0 = flmbd[-1]
    else:
        pi0=1
    if verbose and estimatePi0:
        M.plot(lmbd, pi0lmbd)
        M.plot(x,flmbd)
        M.legend(["pi0(lambda", "fit, len(lmbd)==%i"%len(lmbd)])
        M.title("pi0: %.4f"%pi0)
        M.show()
    args = Numeric.argsort(pValsMA)
    q = Numeric.ones(pValsMA.shape, Numeric.Float)
    q[args[m-1]] = pi0 * pValsMA[args[m-1]]
    for i in range(m-1, 0, -1):
        q[args[i-1]] = min(pi0*m*pValsMA[args[i-1]]/i, q[args[i]])
    MA.put(qVals, putInd, q)
    return qVals



if __name__=="__main__":

    #====================== MULTIVARIATE LINEAR REGRESSION ==========================

    # MLR: compare with orange.LinRegLearner
##    import Dicty.DData, Meda.Preproc, orange
##    DD = Dicty.DData.DData_Nancy()
##    p = DD.getMean2d("pkaC")
##    po = Meda.Preproc.ma2orng(p, Meda.Preproc.getTcDomain(13))
##    po.domain.classVar = po.domain.attributes[-1]
##    lr = orange.LinRegLearner(po)
##    pn = Numeric.array(p)
##    mlr = MultLinReg(pn[:,0:12], pn[:,12])

    # mlr: simple test
##    xy1 = Numeric.array([[1,1.1],[0,-0.01],[-1,-1.1],[2,2.3],[-2,-1.5]], Numeric.Float)
##    ml1 = MultLinReg(xy1[:,:-1],xy1[:,-1])
##    et1 = Meda.Preproc.ma2orng(xy1, Meda.Preproc.getTcDomain(2,False,[]))
##    et1.domain.classVar = et1.domain.attributes[-1]
##    lr1 = orange.LinRegLearner(et1)

    # mlr: test with random numbers
##    import numpy.oldnumeric.random_array as RandomArray
##    xy2 = RandomArray.random((5,4))
##    ml2 = MultLinReg(xy2[:,:-1],xy2[:,-1],1)
##    et2 = Meda.Preproc.ma2orng(xy2, Meda.Preproc.getTcDomain(4,False,[]))
##    et2.domain.classVar = et2.domain.attributes[-1]
##    d3 = orange.Domain(et2.domain.attributes[:-1], et2.domain.attributes[-1])
##    et3 = orange.ExampleTable(d3, et2)
##    lr2 = orange.LinRegLearner(et3)


    #====================== 2 WAY ANOVA ==========================

    def print3d(arr3d):
        for k in range(arr3d.shape[1]):
            print "----",k,"----"
            for j in range(arr3d.shape[2]):
                for i in range(arr3d.shape[0]):
                    print arr3d[i,k,j]

##    import Dicty.DData
##    DN = Dicty.DData.DData_Nancy()
##    rawN = DN.getRaw3dAll()
##    # subset of data
##    rawSub = rawN[0,0:3,0:6]
##    rgiSub = [[0,1,2],[3,4,5]]

##    # compare traditional ANOVA and ANOVA using linear regression on a balanced design
##    as = Anova2way(rawSub, 3)
##    print "\ntrad.anova, interact", as.FAprob, as.FBprob, as.FABprob
##    ar0 = Anova2wayLR(rawSub, rgiSub, 0)
##    print "regr.anova, no inter", ar0.FAprob, ar0.FBprob
##    ar1 = Anova2wayLR(rawSub, rgiSub, 1)
##    print "regr.anova, interact", ar1.FAprob, ar1.FBprob, ar1.FABprob
##    
##    # missing data, no empty cells
##    dSubM1 = MA.array(rawSub)
##    dSubM1[-1,-1] = MA.masked
##    arM1_0 = Anova2wayLR(dSubM1, rgiSub, 0)
##    print "\nregr.anova, no inter", arM1_0.FAprob, arM1_0.FBprob
##    arM1_1 = Anova2wayLR(dSubM1, rgiSub, 1)
##    print "regr.anova, interact", arM1_1.FAprob, arM1_1.FBprob, arM1_1.FABprob
##
##    # missing data and empty cell    
##    dSubM3 = MA.array(rawSub)
##    dSubM3[-1,-3:] = MA.masked
##    arM3_0 = Anova2wayLR(dSubM3, rgiSub, 0)
##    print "\nregr.anova, no inter", arM3_0.FAprob, arM3_0.FBprob
##    arM3_1 = Anova2wayLR(dSubM3, rgiSub, 1)
##    print "regr.anova, interact", arM3_1.FAprob, arM3_1.FBprob, arM3_1.FABprob

    # test LR ANOVA for selection of genes on a whole dataset
    import time

    def anova_2fact_testSpeed(rawAll, rgi):
        """factors: time, strain
        """
        as1 = Numeric.zeros((rawAll.shape[0],2), Numeric.Float)
        print "2 fact start", time.ctime()
        for i in range(rawAll.shape[0]):
            a0 = Anova2wayLR(rawAll[i], rgi, 0)
            as1[i] = [a0.FAprob, a0.FBprob]
        print "2 fact end  ", time.ctime()
        return as1
        
    def anova_fulfact_testSpeed(rawAll, rgi):
        """factors: time, strain & their interaction
        """
        try:
            as1 = Numeric.zeros((rawAll.shape[0],3), Numeric.Float)
            print "2 fact start", time.ctime()
            for i in range(rawAll.shape[0]):
                a0 = Anova2wayLR(rawAll[i], rgi, 1)
                as1[i] = [a0.FAprob, a0.FBprob, a0.FABprob]
            print "2 fact end  ", time.ctime()
        except:
            print "gene", i, time.ctime()
            raise
        return as1

##    asff = anova_fulfact_testSpeed(rawN, DN.replicaGroupInd)  # takes appx. 3 hours!
##    as2f = anova_2fact_testSpeed(rawN, DN.replicaGroupInd)      # takes appx 2 minutes

    #====================== 3 WAY ANOVA ==========================
    
##    import Dicty.DData
##    DN = Dicty.DData.DData_Nancy()
##    rawN = DN.getRaw3dAll()

    # small test: Anova3way
##    rawSub = rawN[0:3, 0:3, 0:3]
##    a3Sub = Anova3way(rawSub, True)
##    print "FprobGene, FprobTime, FprobRep", a3Sub.FAprob, a3Sub.FBprob, a3Sub.FCprob

    # are replicas equal? if FprobRep == 0: different, else equal
##    DN = Dicty.DData.DData_Nancy()
##    DDPJ = Dicty.DData.DData_PJ1()
##    DDMEDSTD = Dicty.DData.DData_MedStd()
##    def anova_testRepEq(DD, addInterAB):
##        print
##        for st in DD.strains:
##            a3 = Anova3way(DD.getRaw3dStrains([st]).filled(0),addInterAB)
##            print "FpG %.10f FpT %.10f FpGxT %.10f FpRep(F) %.20f(%3.5f)\t%12s (%i rep)" % (a3.FAprob, a3.FBprob, a3.FABprob, a3.FCprob, a3.FC, st, len(DD.replicaGroupInd[DD.getIdxStrain(st)]))
##
##    print "------------------ START ---------------"
##    anova_testRepEq(DN, True)
##    anova_testRepEq(DDPJ, True)
##    anova_testRepEq(DDMEDSTD, True)
##    print "------------------ END ---------------"

    #====================== 1 WAY ANOVA LR ==========================

##    a1 = Anova1wayLR(Numeric.ravel(rawSub), [range(0,10), range(10,15), range(15,27)])
##    subT = MA.transpose(rawN[0,:,8:11])

    #====================== 1 WAY REPEATED MEASURES ANOVA LR ==========================

##    a1rm = AnovaRM11LR(subT)
##    a1rmf = AnovaRM11LR(subT.filled(-1))

    #====================== 2 WAY REPEATED MEASURES (on factor A) ANOVA LR ==========================

##    ab = AnovaLRBase()
##    print
##    print Numeric.concatenate((a2rm1.dummyA,a2rm1.dummyB,a2rm1.dummySubjInB),1)
##    print
##
##    import Dicty.DData
##    DN = Dicty.DData.DData_Nancy()
##    rawN = DN.getRaw3dAll()
##    rawSub = rawN[0, 0:5, 0:9]
##    gi = [[0,1,2],[3,4,5],[6,7,8]]

##    a2rm1 = AnovaRM12LR(rawSub,gi,0)
##    print "---no missing---"
##    print a2rm1.FAprob, a2rm1.FBprob#, a2rm1.FABprob
##
##    rawMiss = MA.array(rawSub)
##    rawMiss[0,-1] = MA.masked
##    a2m = AnovaRM12LR(rawMiss,gi,1)
##    print "---with missing---"
##    print a2m.FAprob, a2m.FBprob, a2m.FABprob


    #====================== q-values ==========================

    import matplotlib.pylab as M

    for st, epist in epist_MMAG_nonRM.items():
        for i in range(3):
            print st, i
            ps1 = epist.posthoc_pIntByStrains[:,i].compressed()
            qs1 = qVals(ps, "loess", True)
            ps2 = epist.posthoc_pIntByStrains[:,i].compressed()
            qs2 = qVals(ps2, "loess", True)
            M.title("%s, %i" %(st,i))
            M.subplot(1,2,1)
            M.scatter(ps1, qs1, 0.001)
            M.subplot(1,2,2)
            M.scatter(ps2, qs2, 0.001)
            M.show()

     # increasing pi0lmbd function
##    qqc = qVals(ppc)

     # negative pi0
##    pp = epist_DS_nonRM["pkaRregA"].posthoc_pIntByStrains[:,0].compressed()
##    qq = qVals(pp)
##    M.scatter(pp,qq)
##    M.show()
    
