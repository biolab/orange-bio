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
        
    def _logbin(self, n ,k):
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
        
        
            