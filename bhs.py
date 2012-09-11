from sage.symbolic.function import BuiltinFunction, is_inexact 

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
    
    def _eval_(self, *args):
        if len(args) != self.m + self.n + 2:
            raise RuntimeError, 'No correct input %s' % args

        self.list_a = args[0:self.m]
        self.list_b = args[self.m:self.m+self.n]
        self.q = args[-2]
        self.z = args[-1]

        self.has_arguments = True

        if self.expand:
            return self.z
        return None

    def _evalf_(self, *args, **kwds):
        return 1.0

print "Small test file"
a, b, q, z = var('a b q z')
bhs1 = BasicHypergeometricSeries(1, 1, expand=False)
bhs2 = bhs1(a, b, q, z)
