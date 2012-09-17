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

## Clebsch-Gordan coefficients
class ClebschGordanCoefficients():
    '''Clebsch-Gordan Coefficients for SU_q(2):
    (i)   l_1 - l_2 \leq j \leq l_2 - l_1 \leq l \leq l_1 + l_2
    (ii)  l_2 - l_1 \leq j \leq l_1 - l_2 \leq l \leq l_1 + l_2
    (iii) j \leq l_1 - l_2 \leq -j \leq l \leq l_1 + l_2
    (iv)  -j \leq l_1 - l_2 \leq j \leq l \leq l_1 + l_2'''

    def __init__(self, q, l1, l2, l, j1, j2, j, expand=True):
        self.l1 = l1
        self.l2 = l2
        self.l = l
        self.q = q
        self.j1 = j1
        self.j2 = j2
        self.j = j
        self.expand = expand

    def __repr__(self):
        return 'C_{%s}(%s, %s, %s; %s, %s, %s)' % (self.q, self.l1, self.l2, self.l, \
                                                   self.j1, self.j2, self.j)

    def _latex_(self):
        return 'C^{%s, %s, %s}_{%s; %s %s %s}' % (self.l1, self.l2, self.l, self.q, \
                                                  self.j1, self.j2, self.j)

    def transform1(self):
        '''Transform: (i) <=> (iv) and (ii) <=> (iii)'''
        l1, l2, j1, j2 = self.l1, self.l2, self.j1, self.j2
        
        self.l1 = (l1 + l2 - self.j)/2
        self.l2 = (l1 + l2 + self.j)/2
        self.j1 = (-l1 + l2 + j1 - j2)/2
        self.j2 = (-l1 + l2 - j1 + j2)/2
        self.j = l2 - l1

    def transform2(self):
        '''Transform: (i) <=> (ii) and (iii) <=> (iv)'''
        l1, l2, j1, j2 = self.l1, self.l2, self.j1, self.j2
        
        self.l1 = l2
        self.l2 = l1
        self.j1 = -j2
        self.j2 = -j1
        self.j  = -self.j

    def currentClass(self):
        l1, l2, l, j1, j2, j = self.l1, self.l2, self.l, self.j1, self.j2, self.j

        if abs(j1) > l1 or abs(j2) > l2 or l > l1+l2:
            return 0    # No class

        if l1-l2 <= j and j <= l2-l1 and l2-l1 <= l:
            return 1
        elif l2-l1 <= j and j <= l1-l2 and l1-l2 <= l:
            return 2
        elif j <= l1-l2 and l1-l2 <= -j and -j <= l:
            return 3
        elif -j <= l1-l2 and l1-l2 <= j and j <= l:
            return 4

        return 0

    def transform2one(self):
        classnr = self.currentClass()
        if classnr == 2:
            self.transform2()
        elif classnr == 3:
            self.transform1()
            self.transform2()
        elif classnr == 4:
            self.transform1()
        elif classnr == 0:
            return 0
    
    def evaluate(self):
        if self.currentClass() == 0:
            return 0

        qpoch = qPochhammer(1, expand=self.expand)
        bhs   = BasicHypergeometricSeries(3,2,expand=self.expand)
        l1, l2, l, q, j1, j2, j = self.l1, self.l2, self.l, self.q, \
                                  self.j1, self.j2, self.j

        # line 1 + end of line 3
        sign = (-1)**(2*l1-l2+l+j1) * q**(-2*l1*(2*l1+1)-l2+l+j1-(l2-l1-j)*(l+l1+j-j1))
        qpoch1 = qpoch(q**2, q**2, 2*l1)
        qpoch2 = qpoch(q**2, q**2, -l1+l2-j)
        line1 = sign*qpoch1/qpoch2

        # line 2
        qpoch3 = qpoch(q**(-2*(1+2*l)), q**2, 1)
        qpoch4 = qpoch(q**(-2), q**(-2), l2-j2)
        qpoch5 = qpoch(q**(-2), q**(-2), l2+j2)
        qpoch6 = qpoch(q**(-2), q**(-2), l1+l2-l)
        qpoch7 = qpoch(q**(-2), q**(-2), l1+l2+l+1)
        qpoch8 = qpoch(q**(-2), q**(-2), l1-j1)
        qpoch9 = qpoch(q**(-2), q**(-2), l1+j1)
        line2 = sqrt(qpoch3*qpoch4*qpoch5/(qpoch6*qpoch7*qpoch8*qpoch9))

        # line 3
        qpoch10 = qpoch(q**2, q**2, l-l1+l2)
        qpoch11 = qpoch(q**2, q**2, l-j)
        qpoch12 = qpoch(q**2, q**2, l+l1-l2)
        qpoch13 = qpoch(q**2, q**2, l+j)
        line3 = sqrt(qpoch10*qpoch11/(qpoch12*qpoch13))

        # line 4
        line4 = bhs(q**(-2*(l1-l2+l)), q**(2*(l-l1+l2+1)), q**(-2*(l1+j1)), \
                    q**(2*(l2-l1-j+1)), q**(-4*l1), q**2, q**(-2*(l2+j2)))
        
        
        return line1*line2*line3*line4

# Calculate some simple Clebsch-Gordan coefficients
#q = var('q')
#p = 1
#l = 1
#for i in range(-(l+p), (l+p)+1, 2):
#    for j in range(-(l-p), (l-p)+1, 2):
#        cg = ClebschGordanCoefficients(q, (l+p)/2, (l-p)/2, l, i/2, j/2, (i+j)/2)
#        print cg, '=',
#        cg.transform2one()
#        print cg.evaluate()
