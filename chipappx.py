## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

import math
import numpy.oldnumeric as Numeric, numpy.oldnumeric.mlab as MLab
import numpy.oldnumeric.linear_algebra as LinearAlgebra
import scipy.stats

class ApproxOrthPolyBasis:
    """Approximation of expression profiles with orthogonal polynomials."""

    # for filling missing values in F statistcs (0->p=0, 1e9->p=1)
    _F_fillValue = 1e9

    def __init__(self, points, k=None):
        """sets the points for which the k polynomials of max. degree=k-1 are orthogonal
        default k = number of points
        """
        self.points = points
        if k == None:
            self.k = len(points) #/ 2
        else:
            self.k = k
        self.basis = OrthPolyBasis(self.points, self.k)


    def getAppxCoef(self, arr2d, numCoef=None):
        """returns approx. coeffients given the number of coefficients different from zero"""
        assert len(arr2d.shape) == 2        
        if numCoef == None:
            numCoef = self.k
        else:
            assert 0 <= numCoef <= self.k, "number of coefficients not in range [0," + str(self.k) + "]"
        coef = self.basis.getApproxCoeff(arr2d)
        coef[:,numCoef:] = Numeric.zeros((coef.shape[0], coef.shape[1] - numCoef), Numeric.Float)
        return coef


    def coef_maxCut(self, appxCoef):
        """returns the coefficients different from zero up to the abs. max. coefficient
        where the first coefficient is excluded from finding the max.
        accepts 2d matrix of coefficients where rows represent different curves
        """
        assert len(appxCoef.shape) == 2
        k = Numeric.shape(appxCoef)[1]
        maxInd = Numeric.argmax(Numeric.absolute(appxCoef[:,1:]),1) + 1
        lowDiagOnes = Numeric.fromfunction(lambda i,j: i>=j, (k,k))
        coefSelector = Numeric.take(lowDiagOnes, maxInd, 0)
        return appxCoef*coefSelector


    def getAppxCurve(self, appxCoef, points=None):
        """returns polynomial curve approximation evaluated at given points and given the approx. coefficients
        accepts 2d matrix of coefficients where rows represent different curves
        """
        assert len(appxCoef.shape) == 2
        assert points == None or len(points) > 0, "None or list of points expected"
        return self.basis.evalApproxPoly(appxCoef, points)


    def getAppxCoef3d(self, arr3d, numCoef=None, maxCut=False):
        """returns given number of approx. coefficeints cutted-off at the maximal value (excluding the first one)
        """
        assert len(arr3d.shape) == 3
        if numCoef == None:
            numCoef = self.k
        else:
            assert 0 < numCoef <= self.k, "number of coefficients not in range [1," + str(self.k) + "]"
        coef3d = Numeric.zeros((arr3d.shape[0], self.k, arr3d.shape[2]), Numeric.Float)
        if maxCut:
            for idx2 in range(arr3d.shape[2]):
                coef3d[:,:,idx2] = self.coef_maxCut(self.getAppxCoef(arr3d[:,:,idx2], numCoef=numCoef))
        else:
            for idx2 in range(arr3d.shape[2]):
                coef3d[:,:,idx2] = self.getAppxCoef(arr3d[:,:,idx2], numCoef=numCoef)
        return coef3d


    def getAppxCurve3d(self, appxCoef3d, points=None):
        """returns appx. curves given 3d appx. coefficients in axis 1
        """
        assert len(appxCoef3d.shape) == 3
        assert points == None or len(points) > 0, "None or list of points expected"
        if points == None:
            points = self.points
        curve3d = Numeric.zeros((appxCoef3d.shape[0], len(points), appxCoef3d.shape[2]), Numeric.Float)
        for idx2 in range(appxCoef3d.shape[2]):
            curve3d[:,:,idx2] = self.getAppxCurve(appxCoef3d[:,:,idx2], points)
        return curve3d


    def getAppxCoef2d_pVals(self, arr2d, maxNumCoef):
        """Returns[:,j] the probability that coef. j and all subsequent coef. are equal to 0.
        Reference: Ott, pp.606.
        Use only the first k appx. coef. with pvals below some alpha.
        
        The significance of the coef. estimated by comparing the variance drop of the model2 with the variance of the model1, where:
          model 1: complete model with maxNumCoef coef. different from 0
          model 2: reduced model with coef. in range (0,k) different from 0
        null hypothesis (H0): coef. k,k+1,...,maxNumCoef are equal to 0
            if H0 rejected (pval below some alpha (e.g. 0.05) -> there exist an important coef. among coef. k,k+1,...,maxNumCoef
        repeat the test for k=0,1,...maxNumCoef-1
        """
        assert len(arr2d.shape) == 2, "2d array expected"
        assert 0 < maxNumCoef <= self.k
        coefMax = self.getAppxCoef(arr2d, maxNumCoef)
        curveMax = self.getAppxCurve(coefMax)
        SSE1 = Numeric.add.reduce((arr2d-curveMax)**2,1)
        MSE1 = SSE1 / (arr2d.shape[1]-maxNumCoef)
        #print "SSE1[0], MSE1[0]",SSE1[0], MSE1[0]
        pvals = Numeric.zeros((arr2d.shape[0], maxNumCoef), Numeric.Float)
        for k in range(maxNumCoef):  # test cofInd: [maxNum-1, maxNum-2, ..., minNum]
            #print "Keeping %i coeff" % (k)
            shpk = list(coefMax.shape)
            shpk[1] -= k
            coefk = Numeric.concatenate((coefMax[:,:k], Numeric.zeros((shpk), Numeric.Float)),1)
            curvek = self.getAppxCurve(coefk)
            SSE2 = Numeric.add.reduce((arr2d-curvek)**2,1)
            MSdrop =(SSE2-SSE1) / (maxNumCoef-k)
            F = MSdrop / MSE1
            #2007-10-11: F -> F.filled(???)
            pvals[:,k] = scipy.stats.fprob((maxNumCoef-k), arr2d.shape[1]-maxNumCoef, F.filled(ApproxOrthPolyBasis._F_fillValue))
        return pvals


    def getAppxCoef2d_significant(self, arr2d, maxNumCoef, alpha):
        """Returns[:,j] approx. coeffcients with p-value < alpha; subsequent coef. are equal to 0.
        Reference: Ott, pp.606.
        
        The significance of the coef. estimated by comparing the variance drop of the model2 with the variance of the model1, where:
          model 1: complete model with maxNumCoef coef. different from 0
          model 2: reduced model with coef. in range (0,k) different from 0
        null hypothesis (H0): coef. k,k+1,...,maxNumCoef are equal to 0
            if H0 rejected (pval below some alpha (e.g. 0.05) -> there exist an important coef. among coef. k,k+1,...,maxNumCoef
        repeat the test for k=0,1,...maxNumCoef-1
        """
        assert len(arr2d.shape) == 2, "2d array expected"
        assert 0 < maxNumCoef <= self.k
        coefMax = self.getAppxCoef(arr2d, maxNumCoef)
        curveMax = self.getAppxCurve(coefMax)
        SSE1 = Numeric.add.reduce((arr2d-curveMax)**2,1)
        MSE1 = SSE1 / (arr2d.shape[1]-maxNumCoef)
        #print "SSE1[0], MSE1[0]",SSE1[0], MSE1[0]
        pvals = Numeric.zeros((arr2d.shape[0], maxNumCoef), Numeric.Float)
        for k in range(maxNumCoef):  # test cofInd: [maxNum-1, maxNum-2, ..., minNum]
            #print "Keeping %i coeff" % (k)
            shpk = list(coefMax.shape)
            shpk[1] -= k
            coefk = Numeric.concatenate((coefMax[:,:k], Numeric.zeros((shpk), Numeric.Float)),1)
            curvek = self.getAppxCurve(coefk)
            SSE2 = Numeric.add.reduce((arr2d-curvek)**2,1)
            MSdrop =(SSE2-SSE1) / (maxNumCoef-k)
            F = MSdrop / MSE1
            #2007-10-11: F -> F.filled(???)
            pvals[:,k] = scipy.stats.fprob((maxNumCoef-k), arr2d.shape[1]-maxNumCoef, F.filled(ApproxOrthPolyBasis._F_fillValue))
        pvals = Numeric.where(pvals > alpha, Numeric.resize(Numeric.arange(pvals.shape[1]),pvals.shape), pvals.shape[1])    # MAX where significant, idx where nonsignificant
        firstNonSignIdx = MLab.min(pvals, 1)    # idx of the first non-significant coef.
        coefSign = Numeric.zeros(coefMax.shape, Numeric.Float)
        for idx in range(coefSign.shape[1]):
            coefSign[:,idx] = Numeric.where(idx<firstNonSignIdx, coefMax[:,idx], 0)
        return coefSign
        
    

class OrthPolyBasis:
    """Orthogonal polynomials on a set of distinct points of degree 0 to k-1.
    - best satisfied for x=Numeric.arange(0,4,step) for arbitrary step
    - fragile with increasing k and increasing range of x -> for large i T[i] tends to overflow
    - orthonormal basis not fragile to increasing range of x
    - points mapped to [-1,1]
    TODO: plot() should be smoother !!!
    DONE: go to interval [-1,1], test that k=100 works!!!
    """
    NORM_NONE = 0           # do not normalize the basis polynomials
    NORM_NORM = 1           # generate orthonormal basis: divide base polynomial by sqrt sum squared values at approx. points
    NORM_NORM_T0_1 = 2      # same as NORM_NORM + the first polynomial = 1
    NORM_END1 = 3           # value of all polynomials at the last point equals 1

    def __init__(self, points, k, normalization=NORM_NORM_T0_1, force=False):
        """
        calculate k polynomials of degree 0 to k-1 orthogonal on a set of distinct points
        map points to interval [-1,1]
        INPUT:  points: array of dictinct points where polynomials are orthogonal
                k: number of polynomials of degree 0 to k-1
                force=True creates basis even if orthogonality is not satisfied due to numerical error
        USES:   x: array of points mapped to [-1,1]
                T_: matrix of values of polynomials calculated at x, shape (k,len(x))
                TT_ = T_ * Numeric.transpose(T_)
                TTinv_ = inverse(TT_)
                sc_: scaling factors
                a, b: coefficients for calculating T (2k-4 different from 0, i.e. 6 for k=5)
                n: number of points = len(points)
                normalization = {0|1|2}
        """
        self.k = k  # number of basis polynomials of order 0 to k-1
        self._force = force
        self.points = Numeric.asarray(points, Numeric.Float)
        self.pointsMin = min(points)
        self.pointsMax = max(points)
        # scaling x to [-1,1] results in smaller a and b, T is not affected; overflow is NOT a problem!
        self.xMin = -1
        self.xMax = 1
        self.x = self._map(self.points, self.pointsMin, self.pointsMax, self.xMin, self.xMax)
        # calculate basis polynomials
        self.n = len(points) # the number of approximation points
        t = Numeric.zeros((k,self.n),Numeric.Float)
        a = Numeric.zeros((k,1),Numeric.Float)
        b = Numeric.zeros((k,1),Numeric.Float)
        t[0,:] = Numeric.ones(self.n,Numeric.Float)
        if k > 1: t[1,:] = self.x - sum(self.x)/self.n
        for i in range(1,k-1):
            a[i+1] = Numeric.innerproduct(self.x, t[i,:] * t[i,:]) / Numeric.innerproduct(t[i,:],t[i,:])
            b[i] = Numeric.innerproduct(t[i,:], t[i,:]) / Numeric.innerproduct(t[i-1,:],t[i-1,:])
            t[i+1,:] = (self.x - a[i+1]) * t[i,:] - b[i] * t[i-1,:]
        self.a = a
        self.b = b
        # prepare for approximation
        self._T0 = t
        # orthonormal
        _TT0 = Numeric.matrixmultiply(self._T0, Numeric.transpose(self._T0))
        self.sc1 = Numeric.sqrt(Numeric.reshape(Numeric.diagonal(_TT0),(self.k,1))) # scaling factors = sqrt sum squared self._T0
        self._T1 = self._T0 / self.sc1
        # orthonormal and T[0] == 1
        self.sc2 = Numeric.sqrt(Numeric.reshape(Numeric.diagonal(_TT0),(self.k,1)) / self.n) # scaling factors = sqrt 1/n * sum squared self._T0
        self._T2 = self._T0 / self.sc2
        # T[:,-1] == 1
        self.sc3 = Numeric.take(self._T0, (-1,), 1) # scaling factors = self._T0[:,-1]
        self._T3 = self._T0 / self.sc3
        # set the variables according to the chosen normalization
        self.setNormalization(normalization)


    def getNormalization(self):
        return self._normalization


    def setNormalization(self, normalization):
        if normalization == OrthPolyBasis.NORM_NONE:
            self.T = self._T0
        elif normalization == OrthPolyBasis.NORM_NORM:
            self.T = self._T1
        elif normalization == OrthPolyBasis.NORM_NORM_T0_1:
            self.T = self._T2
        elif normalization == OrthPolyBasis.NORM_END1:
            self.T = self._T3
        else:
            raise "Error: unknown normalization: " + str(normalization)
        self.TT = Numeric.matrixmultiply(self.T, Numeric.transpose(self.T))
        self.TTinv = LinearAlgebra.inverse(self.TT)
        self.TTinvT = Numeric.matrixmultiply(self.TTinv, self.T)
        self.basisCoef = self._getBasisCoef(self.x, self.T)
        self._normalization = normalization
        self._checkOrth(self.T, self.TT, output = self._force)
##        print "self.T", self.T
##        print "self.TT", self.TT
##        print "self.TTinv", self.TTinv
##        print "self.TTinvT", self.TTinvT
        
        

    def _checkOrth(self, T, TT, eps=0.0001, output=False):
        """check if the basis is orthogonal on a set of points x: TT == T*transpose(T) == c*Identity
        INPUT:  T: matrix of values of polynomials calculated at common reference points (x)
                TT = T * transpose(T)
                eps: max numeric error
        """
        TTd0 = (-1.*Numeric.identity(Numeric.shape(TT)[0])+1) * TT  # TTd0 = TT with 0s on the main diagonal
        s = Numeric.sum(Numeric.sum(Numeric.absolute(TTd0)))
        minT = MLab.min(MLab.min(T))
        maxT = MLab.max(MLab.max(T))
        minTTd0 = MLab.min(MLab.min(TTd0))
        maxTTd0 = MLab.max(MLab.max(TTd0))
        if not s < eps:
            out = "NOT ORTHOG, min(T), max(T):\t%f\t%f\n" % (minT, maxT)
            out += "  min(TTd0), max(TTd0), sum-abs-el(TTd0):\t%f\t%f\t%f" % (minTTd0, maxTTd0, s)
            if output:
                print out
                return False
            else:
                raise out
        elif output:
            out = "ORTHOGONAL, min(T), max(T):\t%f\t%f\n" % (minT, maxT)
            out += "  min(TTd0), max(TTd0), sum-abs-el(TTd0):\t%f\t%f\t%f" % (minTTd0, maxTTd0, s)
            print out
        return True


    def zipab(self):
        return zip(Numeric.reshape(self.a,(len(self.a),)).tolist(),Numeric.reshape(self.b, (len(self.b),)).tolist())


    def _map(self, x, a,b, m,M):
        """maps x from the interval [a,b] to the interval [m,M]"""
        x = Numeric.asarray(x, Numeric.Float)
        return (M-m)*(x-a)/(b-a) + m


    #--------APPROXIMATION-------------------------


    def getApproxCoeff(self, curves):
        """curves: 1d or 2d array where each curve in separate line, e.g. curves[curveIdx,timePoints]"""
##        print "curves"
##        print curves
##        print
##        print "Numeric.transpose(self.TTinvT)"
##        print Numeric.transpose(self.TTinvT)
##        print
##        print "Numeric.matrixmultiply(curves, Numeric.transpose(self.TTinvT))"
##        print Numeric.matrixmultiply(curves, Numeric.transpose(self.TTinvT))
        return Numeric.matrixmultiply(curves, Numeric.transpose(self.TTinvT))


    #--------EVALUATION OF POLYNOMIAL-------------------------

    def _getBasisCoef(self, x, T):
        """returns coefficients of basis polynomials given T, i.e. their values at x: (poly0, xi)
        where polynomials in rows, coefficients [a0,a1,...ak] where k=len(x)-1
        similar goes for T
        """
        numPoints = T.shape[1]  # number of points that polynomials are calculated == len(x)
        assert len(x) == numPoints, "len(x)=" + str(x) + ", T.shape[1]=" + str(numPoints) + " do not match"
        numLinSys = T.shape[0]  # number of polynomials
        lin_sys_b = Numeric.transpose(T)
        lin_sys_a = Numeric.ones((numPoints, numPoints), Numeric.Float)
        lin_sys_a[:,1] = x
        for idx in range(2,numPoints):
            lin_sys_a[:,idx] = lin_sys_a[:,idx-1] * x
        return Numeric.transpose(LinearAlgebra.solve_linear_equations(lin_sys_a, lin_sys_b))


    def _evalPolyHorner(self, polyCoef, x):
        """returns (#poly, #x) polynomials evaluated at x where:
            - result: polynomials in rows, values at x in columns: [f(x0)...f(xn)]
            - polyCoef: polynomials in rows, coefficients in columns, sorted by increasing power of x: [c0, c1, ...]
        uses Horner's rule for polynomial evaluation
        """
        x = Numeric.asarray(x, Numeric.Float)
        one_poly = len(polyCoef.shape) == 1
        if one_poly:
            polyCoef = polyCoef[Numeric.NewAxis,:] # [polyIdx, coefIdx]
        val = Numeric.zeros((polyCoef.shape[0],len(x)), Numeric.Float)  # [polyIdx, xIdx]
        for idx in range(polyCoef.shape[1]-1,0,-1):
            val = (val + polyCoef[:,idx:idx+1]) * x
        val += polyCoef[:,0:1]
        if one_poly:
            return val[0,:]
        else:
            return val


    def evalApproxPoly(self, appxCoef, points=None):
        """returns (#curve, #points) an approx. polynomials calculated at points given approx. coeff in rows:
            - appxCoef: curves for approximation in rows, appxCoef in columns [B0, B1, ...]
            - points relative to self.points
        TODO: evaluate on a matrix of appxCoef
        """
        if points == None:
            return Numeric.matrixmultiply(appxCoef, self.T)
        appxCoef = Numeric.asarray(appxCoef, Numeric.Float)
        one_curve = len(appxCoef.shape) == 1
        if one_curve:
            appxCoef = appxCoef[Numeric.NewAxis,:]  # [curveIdx, coefIdx]
        mappedPoints = self._map(points, self.pointsMin, self.pointsMax, self.xMin, self.xMax)
        # eval basis polynomials on mapped points
        basisEval = self._evalPolyHorner(self.basisCoef, mappedPoints) #[basisIdx == coefIdx, pointIdx]
        result = Numeric.matrixmultiply(appxCoef, basisEval)
        if one_curve:
            return result[0,:]
        else:
            return result


class TrigonomerticBasis:
    """Approximation of expression profiles with trigonometric functions orthogonal polynomials."""

    def __init__(self, numPoints, k):
        """numPoints: number of approximation points; k: number of basis functions [2,...,numPoints]"""
        self.numPoints = numPoints
        self.k = k
##        assert k > 1, "Error TrigonomerticBasis: k <= 1"
        assert k <= numPoints, "Error TrigonomerticBasis: k > numPoints"
        # evaluate trigonometric basis functions on the given number of points from [-pi,pi]
        self.x = Numeric.arange(-1*math.pi, math.pi+0.0000001, 2*math.pi/(numPoints-1))
        self.y = Numeric.ones((k, numPoints), Numeric.Float)
        for kk in range(1, k, 2):
##            print "kk, cos %ix" % ((kk+1)/2.)
            self.y[kk] = MLab.cos(self.x*(kk+1)/2) 
        for kk in range(2, k, 2):
##            print "kk, sin %ix" % (kk/2.)
            self.y[kk] = MLab.sin(self.x*kk/2)
        # approx. matrix
        self.Ainv = LinearAlgebra.inverse(Numeric.matrixmultiply(self.y, Numeric.transpose(self.y)))
        self.yyTinvy = Numeric.matrixmultiply(LinearAlgebra.inverse(Numeric.matrixmultiply(self.y, Numeric.transpose(self.y))), self.y)

    def getAppxCoef(self, curves):
        return Numeric.matrixmultiply(curves, Numeric.transpose(self.yyTinvy))

    def getAppxCurve(self, appxCoef):
        return Numeric.matrixmultiply(appxCoef, self.y)

