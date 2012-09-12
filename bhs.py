from sage.symbolic.function import BuiltinFunction, is_inexact 

class qPochhammer(BuiltinFunction):
    
    def __init__(self, m=1, expand=False):
        self.m = m
        self.expand = expand
        self.has_arguments = False

        BuiltinFunction.__init__(self, 'qPochhammer(%i)' % self.m, nargs=(m+2))

    def _latex_(self):
        if self.has_arguments:
            return '(%s; %s)_{%s}' % (self.a, self.q, self.n)
        return 'qPochhammer'

    def _qPochhammer1(self, a, q, n):
        if n == 0:
            return 1
        return (1 - a*q**(n-1))*self._qPochhammer1(a, q, n-1)

    def _expand(self):
        return prod([self._qPochhammer1(elm, self.q, self.n) \
                     for elm in self.a])

    def _eval_(self, *args):
        if len(args) != self.m + 2:
            raise RuntimeError, args
        
        if self.expand:        
            return self.evaluate(args[0:self.m], args[-2], args[-1])
        
        return None

    def _evalf_(self, *args, **kwds):
        if type(self.a) == tuple:
            return self._expand()
        return self._qPochhammer1(self.a, self.q, self.n)
        
    def evaluate(self, a, q, n):
        self.a = a
        self.q = q
        self.n = n

        self.has_arguments = True

        if type(self.a) == tuple:
            return self._expand()
        return self._qPochhammer1(self.a, self.q, self.n)


class BasicHypergeometricSeries(BuiltinFunction):
    
    def __init__(self, m=2, n=1, expand=False):
        self.m = m
        self.n = n
        self.expand = expand
        self.has_arguments = False

        BuiltinFunction.__init__(self, 'bhs(%i,%i)' % (m, n), nargs=(m+n+2), 
                                 latex_name='{}_{%i}\\phi_{%i}' % (m, n))

    def _latex_(self):
        if self.has_arguments:
            pass
        return '{}_{%i}\\phi_{%i}' % (self.m, self.n)
    
    def expand_bhs(self, k=0):
        coef = self._expand_coefficient(k)
        if coef != 0:
            return coef*self.z**k + self._expand(k+1)
        return 0

    def expand_coefficient(self, k):
        nom = qPochhammer(len(self.list_a)).evaluate(self.list_a, self.q, k)
        if nom == 0:
            return 0

        denom = qPochhammer(len(self.list_b)).evaluate(self.list_b, self.q, k)
        
        coef = nom / denom * \
            ((-1)**k) * q**(binomial(n, 2))**(len(self.list_b) - len(self.list_a))

        return coef

    def _eval_(self, *args):
        if len(args) != self.m + self.n + 2:
            raise RuntimeError, 'No correct input %s' % args

        self.list_a = args[0:self.m]
        self.list_b = args[self.m:self.m+self.n]
        self.q = args[-2]
        self.z = args[-1]

        self.has_arguments = True

        if self.expand:
            return self._expand()
        return None

    def _evalf_(self, *args, **kwds):
        return self._expand()

print "Small test file"
a, b, q, z = var('a b q z')
bhs1 = BasicHypergeometricSeries(2, 1, expand=False)
bhs2 = bhs1(q**(-2), a, b, q, z)
