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

        sage: a, b, c, q = var('a b c q')
        sage: qshift3 = qPochhammer(m=3, expand_exp=True)
        sage: qshift3(a, b, c, q, 2)
        (c - 1)*(b - 1)*(a - 1)*(c*q - 1)*(b*q - 1)*(a*q - 1)

    """

    def __init__(self, m=1, expand_exp=False):
        self.m = m
        self.expand_exp = expand_exp

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

        # Implementation for infinity oo
        if args[-1] == oo:
            return None # no expantion
        
        if self.expand_exp:
            return self.evaluate(args[:self.m], args[-2], args[-1])
        
        return None

    def _evalf_(self, *args, **kwds):
        if type(self.a) == tuple:
            return self._expand()
        return self._qPochhammer1(self.a, self.q, self.n)
        
    def evaluate(self, a, q, n):
        
        if type(self.a) == tuple:
            return self._expand()
        return self._qPochhammer1(a, q, n)


class BasicHypergeometricSeries(BuiltinFunction):
    r"""
    BasicHypergeometricSeries computes the Basic Hypergeometric Series

    INPUT:
        - ``m`` - Integer, 2 by default, gives the number of elements in the
                  nominator.
        - ``n`` - Integer, 1 by default, gives the number of elements in the
                  denominator,
        - ``expand_exp`` - Boolean, False by default, if expression should be
                           expanded.

    OUTPUT:
        - Builtinfunction which can compute certain basic hypergeometric 
          series.

    EXAMPLES:

    Most simple example::

        sage: a, b, c, q, z = var('a b c q z')
        sage: phi = BasicHypergeometricSeries()
        sage: phi(a, b, c, q, z)
        bhs(2,1)(a, b, c, q, z)

    Expand expression::

        sage: phi = BasicHypergeometricSeries(expand_exp=True)
        sage: phi(q**(-2), a, b, q, z)
        (1/q^2 - 1)*(1/q - 1)*(a - 1)*(a*q - 1)*z^2/((b - 1)*(b*q - 1)*q) 
            + (1/q^2 - 1)*(a - 1)*z/(b - 1) + 1

    Use different amount of nominator/denominator::

        sage: phi = BasicHypergeometricSeries(m=4, n=2, expand_exp=True)
        sage: phi = BasicHypergeometricSeries(m=4, n=2, expand_exp=True)
        sage: phi(q**(-2), a1, a2, a3, b1, b2, q, z)
        (1/q^2 - 1)*(1/q - 1)*(a3 - 1)*(a2 - 1)*(a1 - 1)*(a3*q - 1)*(a2*q - 1)
            *(a1*q - 1)*z^2/((b2 - 1)*(b1 - 1)*(b2*q - 1)*(b1*q - 1)*q^2) 
            + (1/q^2 - 1)*(a3 - 1)*(a2 - 1)*(a1 - 1)*z/((b2 - 1)*(b1 - 1)) + 1

    """
    
    def __init__(self, m=2, n=1, expand_exp=False):
        self.m = m
        self.n = n
        self.expand = expand_exp
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


