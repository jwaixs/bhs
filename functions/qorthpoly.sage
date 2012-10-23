load('bhs.sage')

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


class qKrawtchoukPolynomials(BuiltinFunction):

    def __init__(self, expand_exp=True):
        self.expand = expand_exp

        BuiltinFunction.__init__(self, 'qkp', nargs=5, latex_name="p")

    def _eval_(self, *args):
        if len(args) != 5:
            raise RuntimeError, args
        
        (n, x, p, N, q) = args
        if n < 0 or n > N:
            raise RuntimeError, 'n = %i should be between 0 and N = %i' % (n, N)
        bhs = BasicHypergeometricSeries(m=3,n=2,expand_exp=self.expand)
        if self.expand:
            return bhs(q**(-n), q**(-x), -p*q**n, q**(-N), 0, q, q)
  
        return None

    def general(self, n, N):
        x, p, q = var('x p q')
        return self._eval_(n, x, p, N, q)

class DualqKrawchoukPolynomials(BuiltinFunction):

    def __init__(self, expand_exp=True):
        self.expand = expand_exp

        BuiltinFunction.__init__(self, 'dqkp', nargs=5, latex_name='dqkp')

    def _eval_(self, *args):
        if len(args) != 5:
            raise RuntimeError, args

        (n, x, c, N, q) = args
        bhs = BasicHypergeometricSeries(m=3, n=2, expand_exp=self.expand)
        if self.expand:
            return bhs(q**(-n), q**(-x), c*q**(x-N), q**(-N), 0, q, q)

        return None

    def general(self, n, N):
        x, c, q = var('x c q')
        return self._eval_(n, x, c, N, q)

# At the moment (10 October 2012) these are unknown polynomials. I'll call them
# Twisted Quantum q-Krawtchouk polynomials.

def twisted_quantum_q_krawtchouk_polynomials(N, q):
    polynomials = [1, x/(1-q**N)]
    for n in range(2, N+1):
        polynomial = x/(1-q**(N-n+1))*polynomials[n-1] - \
                     (1-q**(n-1))/(1-q**(N-n+1))*polynomials[n-2]
        polynomials.append(polynomial)
    return polynomials

#class TwistedQuantumqKrawtchoukPolynomials():
#
#    def __init__(self, N):
#        self.N = N
#        self.evaluate_monic_polynomials(self.N)
#        self.evaluate_easy_polynomials(self.N)
#        
#    def evaluate_monic_polynomials(self, N):
#        q, x = var('q x')
#        result = [1, x]
#
#        for r in range(2,N+1):
#            pr = x*result[r-1] - (1-q**(r-1))*(1-q**(N-r+2))*result[r-2]
#            result.append(pr)
#
#        self.monic_polynomials = result
#
#    def evaluate_easy_polynomials(self, N):
#        q, x = var('q x')
#        result = [1, x/(1-q**N)]
#
#        for r in range(2, N+1):
#            pr = x*result[r-1]/(1 - q**(N-r+1)) \
#                - (1-q**(r-1))/(1 - q**(N-r+1))*result[r-2]
#            result.append(pr)
#
#        self.easy_polynomials = result
#
