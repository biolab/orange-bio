import math

class LogBin(object):
    _max = 2
    _lookup = [0, 1]
    def __init__(self, max=1000):
        self._extend(max)

    def  _extend(self, max):
        lookup = [self._lookup[-1] + math.log(self._max)] * (max - self._max)
        for i in xrange(1, max - self._max):
            lookup[i] = lookup[i - 1] + math.log(self._max + i)
        LogBin._lookup.extend(lookup)
        LogBin._max = max
        
    def _logbin(self, n, k):
        if n >= self._max:
            self._extend(n + 100)
        return self._lookup[n] - self._lookup[n - k] - self._lookup[k]

class Binomial(LogBin):
##    def __init__(self, max=1000):
##        self.max = 2
##        self._lookup = [0, 1]
##        self._extend(max)
##
##    def  _extend(self, max):
##        lookup = [self._lookup[-1] + math.log(self.max)] * (max - self.max)
##        for i in xrange(1, max- self.max):
##            lookup[i] = lookup[i - 1] + math.log(self.max + i)
##        self._lookup.extend(lookup)
##        self.max = max
##        
##    def _logbin(self, n ,r):
##        if n >= self.max:
##            self._extend(n+100)
##        return self._lookup[n] - self._lookup[n-r] - self._lookup[r]    

    def __call__(self, k, N, m, n):
        p = 1.0 * m / N
        if p == 0.0:
            if k == 0:
                return 1.0
            else:
                return 0.0
        elif p == 1.0:
            if n == k:
                return 1.0
            else:
                return 0.0
        return math.exp(self._logbin(n, k) + k * math.log(p) + (n + k) * math.log(1.0 - p))

    def p_value(self, k, N, m, n):
        subtract = False #n - k + 1 > k
        result = sum(self.__call__(i, N, m, n) for i in (range(k) if subtract else range(k, n+1)))
        return 1.0 - result if subtract else result

class Hypergeometric(LogBin):
    def __call__(self, k, N, m, n):
        return math.exp(self._logbin(m, k) + self._logbin(N - m, n - k) - self._logbin(N, n))

    def p_value(self, k, N, m, n):
        subtract = False #n - k + 1 > k
        result = sum(math.exp(self._logbin(m, i) + self._logbin(N - m, n - i)) for i in (range(k) if subtract else range(k, n+1)))
        result /= math.exp(self._logbin(N, n))
        return 1.0 - result if subtract else result


## to speed-up FDR, calculate ahead sum([1/i for i in range(1, m+1)]), for m in [1,100000]. For higher values of m use an approximation, with error less or equal to 4.99999157277e-006. (sum([1/i for i in range(1, m+1)])  ~ log(m) + 0.5772..., 0.5572 is an Euler-Mascheroni constant) 
c = [1.0]
for m in range(2, 100000):
    c.append( c[-1] + 1.0/m)

def FDR(p_values, q=0.05, dependent=False, m=None):
    if not(m):
        m = len(p_values)
    if m == 0:
        return []

    if dependent: # correct q for dependent tests
        k = c[m-1] if m <= len(c) else math.log(m) + 0.57721566490153286060651209008240243104215933593992
        q = q/k

    return [p for (i, p) in enumerate(p_values) if p <= (i+1.0)*q/m]

def Bonferroni(p_values, q=0.05, m=None):
    if not(m):
        m = len(p_values)
    if m == 0:
        return []
    q = q/float(m)
    return [p for p in p_values if p <= q]
