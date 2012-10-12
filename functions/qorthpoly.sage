class littleqJacobiPolynomials(BuiltinFunction):

    def __init__(self, expand_exp=True):
        self.expand = expand_exp
        self.has_arguments = False

        BuiltinFunction.__init__(self, 'lqjp', nargs=5, latex_name="p")

    def _eval_(self, *args):
        if len(args) != 5:
            raise RuntimeError, args
        
        self.x = args[0]
        self.a = args[1]
        self.b = args[2]
        self.q = args[3]
        self.n = args[4]

        self.has_arguments = True

        self.bhs = BasicHypergeometricSeries(m=2,n=1,expand_exp=self.expand)
        if self.expand:
            return self.bhs(self.q**(-self.n), self.a*self.b*self.q**(self.n+1), \
                            self.a*self.q, self.q, self.q*self.x)
  
        return None

    def general(self, n):
        x, a, b, q = var('x a b q')
        return self._eval_(x, a, b, q, n)


# At the moment (10 October 2012) these are unknown polynomials. I'll call them
# Twisted Quantum q-Krawtchouk polynomials.

class TwistedQuantumqKrawtchoukPolynomials():

    def __init__(self, N):
        self.N = N
        self.evaluate_polynomials(self.N)
        
    def evaluate_polynomials(self, N):
        q, x = var('q x')
        result = [1, x]

        for r in range(2,N+1):
            pr = x*result[r-1] - (1-q**(r-1))*(1-q**(N-r+2))*result[r-2]
            result.append(pr)

        self.polynomials = result

