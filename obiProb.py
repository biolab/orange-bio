import math

class Binomial(object):
    def __init__(self, max=1000):
        self.max = 2
        self._lookup = [0, 1]
        self._extend(max)

    def  _extend(self, max):
        lookup = [self._lookup[-1] + math.log(self.max)] * (max - self.max)
        for i in xrange(1, max- self.max):
            lookup[i] = lookup[i - 1] + math.log(self.max + i)
        self._lookup.extend(lookup)
        self.max = max
        
    def _logbin(self, n ,r):
        if n >= self.max:
            self._extend(n+100)
        return self._lookup[n] - self._lookup[n-r] - self._lookup[r]    

    def __call__(self, n, r, p):
        if p==0.0:
            if r==0:
                return 0.0
            else:
                return 1.0
        elif p==1.0:
            if n==r:
                return 0.0
            else:
                return 1.0
        return math.exp(self._logbin(n, r) + r*math.log(p) + (n + r)*math.log(1.0-p))

    def PValue(self, from_, to, p):
        return reduce(lambda sum, i: sum + self.__call__(to, i, p), range(from_, to+1), 0.0)
            