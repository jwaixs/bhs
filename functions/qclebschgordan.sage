load('bhs.sage')

class QuantumClebschGordanCoefficients():
    '''Clebsch-Gordan Coefficients for SU_q(2):
    (i)   l_1 - l_2 \leq j \leq l_2 - l_1 \leq l \leq l_1 + l_2
    (ii)  l_2 - l_1 \leq j \leq l_1 - l_2 \leq l \leq l_1 + l_2
    (iii) j \leq l_1 - l_2 \leq -j \leq l \leq l_1 + l_2
    (iv)  -j \leq l_1 - l_2 \leq j \leq l \leq l_1 + l_2'''

    def __init__(self, q, l1, l2, l, j1, j2, j, expand_exp=True):
        self.l1 = l1
        self.l2 = l2
        self.l = l
        self.q = q
        self.j1 = j1
        self.j2 = j2
        self.j = j
        self.expand = expand_exp

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

        if abs(j1) > l1 or abs(j2) > l2 or l > l1+l2 or j1+j2 != j:
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

        qpoch = qPochhammer(1, expand_exp=True)
        bhs   = BasicHypergeometricSeries(3,2,expand_exp=True)
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

def quantum_clebsch_gordan(q, l1, l2, l, j1, j2, j, algorithm='sage'):
    cg = QuantumClebschGordanCoefficients(q, l1, l2, l, j1, j2, j)
    cg.transform2one()    

    evaluation = cg.evaluate()

    if type(evaluation) != Expression:
        return evaluation

    if algorithm == 'mathematica':
        return mathematica(evaluation).Simplify().sage()
    return evaluation.simplify()
