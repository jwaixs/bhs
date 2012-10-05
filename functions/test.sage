    from sage.symbolic.function import BuiltinFunction

    class qPochhammer(BuiltinFunction):
        
        def __init__(self, m=1, expand_exp=False):
            self.m = m
            self.expand_exp = expand_exp
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
            
            if self.expand_exp:
                return self.evaluate(args[:self.m], args[-2], args[-1])
            
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
