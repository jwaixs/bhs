r"""
Basic Hypergeometric Series

Collection of the most basic operation in Basic Hypergeometric Series. It is
designed to be able to nummerical calculations and algebraic calculations. The
algorithms are based on the computations in [GR04].

REFERENCES:

.. [GR04] G. Gasper, M. Rahman, 'Basic Hypergeometric Series', second edition,
          Encyclopedia of Mathematics and Its Applications 96, Cambridge
          University Press, 2004.

AUTHORS:

- Noud Aldenhoven (2012-10-05): Initial version for Sage
"""

from sage.symbolic.function import BuiltinFunction, is_inexact 
from sage.rings.integer import Integer
from sage.symbolic.expression import Expression

class qPochhammer(BuiltinFunction):
    r"""
    qPochhammer computes the q-shifted factorial.

    INPUT:

    - ``m`` - Integer, 1 by default, number of q-shifted factorials to compute.

    - ``expand_exp`` - Boolean, false by default, if expression should be 
                       expanded.

    OUTPUT:

    - BuiltinFunction which can compute q-shifted factorials.

    EXAMPLES: 

    Most simple example::

        sage: q, a = var('q a')
        sage: qshift = qPochhammer()
        sage: qshift
        qPochhammer(1)
        sage: qshift(a, q, 2)
        qPochhammer(1)(a, q, 2)
 
    Expand expression immediately::

        sage: qshift = qPochhammer(expand_exp=True)
        sage: qshift(a, q, 2)
        (a - 1)*(a*q - 1)

    Expand expression for more q-shifted factorials::

    """

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


class BasicHypergeometricSeries(BuiltinFunction):
    
    def __init__(self, m=2, n=1, expand=False):
        self.m = m
        self.n = n
        self.expand = expand
        self.has_arguments = False

        BuiltinFunction.__init__(self, 'bhs(%i,%i)' % (m, n), nargs=(m+n+2), 
                                 latex_name='{}_{%i}\\phi_{%i}' % (m, n))

    def _latex_(self):
        return '{}_{%i}\\phi_{%i}' % (self.m, self.n)

    def expand_bhs(self, k=0):
        coef = self.expand_coefficient(k)
        if coef != 0:
            return coef*self.z**k + self.expand_bhs(k+1)
        return 0

    def expand_coefficient(self, k):
        nom = qPochhammer(len(self.list_a)).evaluate(self.list_a, self.q, k)
        if nom == 0:
            return 0

        denom = qPochhammer(len(self.list_b)).evaluate(self.list_b, self.q, k)
     
        coef = nom / denom * \
            ((-1)**k * q**(binomial(k, 2)))**(len(self.list_b) - len(self.list_a))

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
            return self.expand_bhs()
        return None

    def _evalf_(self, *args, **kwds):
        return self.expand_bhs()


## q-Orthogonal polynomials
class littleqJacobiPolynomials(BuiltinFunction):

    def __init__(self, expand=True):
        self.expand = expand
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

        self.bhs = BasicHypergeometricSeries(2,1,expand=self.expand)
        if self.expand:
            return self.bhs(self.q**(-self.n), self.a*self.b*self.q**(self.n+1), \
                            self.a*self.q, self.q, self.q*self.x)
  
        return None

    def general(self, n):
        x, a, b, q = var('x a b q')
        return self._eval_(x, a, b, q, n)

