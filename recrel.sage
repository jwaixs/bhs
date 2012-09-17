from itertools import product

load('bhs.sage')

def recrel_a(q, l1, l2, m1, m2, l):
    total_sum = 0
    for i1, i2, j1, j2, n1, n2 in product(range(2*l1+1), range(2*l2+1),\
        range(2), range(2), range(2*m1+1), range(2*m2+1)):
        cg1 = ClebschGordanCoefficients(q, l1, l2, l, i1-l1, i2-l2, l)
        cg2 = ClebschGordanCoefficients(q, 1/2, 1/2, 0, j1-1/2, j2-1/2, 0)
        cg3 = ClebschGordanCoefficients(q, l1, 1/2, m1, i1-l1, j1-1/2, n1-m1)
        cg4 = ClebschGordanCoefficients(q, l2, 1/2, m2, i2-l2, j2-1/2, n2-m2)
        cg5 = ClebschGordanCoefficients(q, m1, m2, l, n1-m1, n2-m2, l)

        cg_list = [cg1, cg2, cg3, cg4, cg5]
        map(lambda elm : elm.transform2one(), cg_list)
        total_sum += prod(map(lambda elm : elm.evaluate(), cg_list))
    
    return total_sum

q = var('q')
print recrel_a(q, 1/2, 1/2, 0, 0, 0).simplify()
print recrel_a(q, 1/2, 1/2, 1, 1, 0).simplify()
