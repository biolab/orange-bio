import math    

class LogBin(object):
    _max = 2
    _lookup = [0.0, 0.0]
    _max_factorial = 1
    def __init__(self, max=1000):
        self._extend(max)

    @classmethod
    def  _extend(cls, max):
        for i in xrange(cls._max, max):
            cls._max_factorial *= i
            cls._lookup.append(math.log(cls._max_factorial))
        cls._max = max

    def _logbin(self, n, k):
        if n >= self._max:
            self._extend(n + 100)
        if k < n:
            return self._lookup[n] - self._lookup[n - k] - self._lookup[k]
        else:
            return 0.0

class Binomial(LogBin):

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
        try:
            return math.exp(self._logbin(n, k) + k * math.log(p) + (n - k) * math.log(1.0 - p))
        except (OverflowError, ValueError), er:
            print k, N, m, n
            raise
##        return math.exp(self._logbin(n, k) + math.log((p**k) * (1.0 - p)**(n - k)))

    def p_value(self, k, N, m, n):
        subtract = n - k + 1 > k
        result = sum([self.__call__(i, N, m, n) for i in (range(k) if subtract else range(k, n+1))])
        return 1.0 - result if subtract else result

class Hypergeometric(LogBin):

    def __call__(self, k, N, m, n):
        try:
            return math.exp(self._logbin(m, k) + self._logbin(N - m, n - k) - self._logbin(N, n))
        except (OverflowError, ValueError), er:
            print k, N, m, n
            raise

    def p_value(self, k, N, m, n):
        subtract = n - k + 1 > k
##        result = sum([math.exp(self._logbin(m, i) + self._logbin(N - m, n - i)) for i in (range(k) if subtract else range(k, n+1))])
        result = sum([self.__call__(i, N, m, n) for i in (range(k) if subtract else range(k, n+1))])
##        result /= math.exp(self._logbin(N, n))
        return 1.0 - result if subtract else result


## to speed-up FDR, calculate ahead sum([1/i for i in range(1, m+1)]), for m in [1,100000]. For higher values of m use an approximation, with error less or equal to 4.99999157277e-006. (sum([1/i for i in range(1, m+1)])  ~ log(m) + 0.5772..., 0.5572 is an Euler-Mascheroni constant) 
c = [1.0]
for m in range(2, 100000):
    c.append( c[-1] + 1.0/m)

def FDR(p_values, dependent=False, m=None):
    if not m:
        m = len(p_values)
    if m == 0:
        return []

    if dependent: # correct q for dependent tests
        k = c[m-1] if m <= len(c) else math.log(m) + 0.57721566490153286060651209008240243104215933593992
        m = m * k

    tmp_fdrs = [p*m/(i+1.0) for (i, p) in enumerate(p_values)]
    fdrs = []
    cmin = tmp_fdrs[-1]
    for f in reversed(tmp_fdrs):
        cmin = min(f, cmin)
        fdrs.append( cmin)
    fdrs.reverse()
    return fdrs

def Bonferroni(p_values, m=None):
    if not m:
        m = len(p_values)
    if m == 0:
        return []
    m = float(m)
    return [p/m for p in p_values]
