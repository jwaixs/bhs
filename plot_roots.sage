# Twisted Dual q-Krawchouk polynomials
#print "Plot zeros for twisted dual q-Krawchouk polynomials"
#roots = load('data/roots')
#
#q_values = [1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8]
max_N = 30
#
#for newq in q_values:
#    points = []
#    for i in range(1, max_N):
#        zeros = roots[newq][i]
#        points += [point(zip(zeros, len(zeros)*[i/max_N]), rgbcolor=hue(i/max_N))]
#    sum(points).save('plots/roots_%s.png' % RR(newq))


# Ordinary Chebychev polynomials of first and second kind
print "Plot zeros for Chebyshev polynomials of the first kind."
x = PolynomialRing(RealField(500), 'x').gen()
points = []
for i in range(1, max_N):
    polynomial = chebyshev_T(i, x)
    zeros = [ zero for (zero, mult) in polynomial.roots() ]
    points += [ point(zip(zeros, len(zeros)*[i/max_N]), rgbcolor=hue(i/max_N)) ]
sum(points).save('plots/chebychev_T.png')

print "Plot zeros for Chebyshev polynomials of the second kind."
x = PolynomialRing(RealField(500), 'x').gen()
points = []
for i in range(1, max_N):
    polynomial = chebyshev_U(i, x)
    zeros = [ zero for (zero, mult) in polynomial.roots() ]
    points += [ point(zip(zeros, len(zeros)*[i/max_N]), rgbcolor=hue(i/max_N)) ]
sum(points).save('plots/chebychev_U.png')
