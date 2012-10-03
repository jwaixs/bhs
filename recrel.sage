from itertools import product

load('bhs.sage')

def recrel_a(q, l1, l2, m1, m2, l):
    total_sum = 0
    for i1, i2, j1, j2, n1, n2 in product(range(2*l1+1), range(2*l2+1),\
        range(2), range(2), range(2*m1+1), range(2*m2+1)):

        if i1-l1+i2-l2 != l or j1+j2-1 != 0 or i1-l1+j1-1/2 != n1-m1 \
                            or i2-l1+j2-1/2 != n2-m2 or n1-m1+n2-m2 != l:
            continue

        cg1 = cgc(q, l1, l2, l, i1-l1, i2-l2, l)
        if sign(cg1.substitute(q=1/2)) != sign(clebsch_gordan(l1, l2, l, i1-l1, i2-l2, l)):
            cg1 = -cg1

        cg2 = cgc(q, 1/2, 1/2, 0, j1-1/2, j2-1/2, 0)
        if sign(cg2.substitute(q=1/2)) != sign(clebsch_gordan(1/2, 1/2, 0, j1-1/2, j2-1/2, 0)):
            cg2 = -cg2

        cg3 = cgc(q, l1, 1/2, m1, i1-l1, j1-1/2, n1-m1)
        if sign(cg3.substitute(q=1/2)) != sign(clebsch_gordan(l1, 1/2, m1, i1-l1, j1-1/2, n1-m1)):
            cg3 = -cg3

        cg4 = cgc(q, l2, 1/2, m2, i2-l2, j2-1/2, n2-m2)
        if sign(cg4.substitute(q=1/2)) != sign(clebsch_gordan(l1, 1/2, m2, i2-l2, j2-1/2, n2-m2)):
            cg4 = -cg4

        cg5 = cgc(q, m1, m2, l, n1-m1, n2-m2, l)
        if sign(cg5.substitute(q=1/2)) != sign(clebsch_gordan(m1, m2, l, n1-m1, n2-m2, l)):
            cg5 = -cg5

        cg_list = [cg1, cg2, cg3, cg4, cg5]
        total_sum += prod(cg_list)

    return total_sum

