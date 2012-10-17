from itertools import product

load('functions/bhs.sage')

q = var('q')
x = var('x')

def recrel_a(q, l1, l2, m1, m2, l):
    r"""
    recrel_a computes the coefficients of Proposition 3.1.

    INPUT:

        ``q``  - where q can be symbolic or 0 < q < 1.
        ``l1`` - rational number where l1 is positive and 2*l1 is an integer.
        ``l2`` - rational number where l2 is positive and 2*l2 is an integer.
        ``m1`` - rational number where -l1 <= m1 <= l1 and 2*m1 is an integer.
        ``m2`` - rational number where -l2 <= m2 <= l2 and 2*m2 is an integer.
        ``l``  - rational number where |l1 - l2| <= l <= l1 + l2 and 2*l is an
                 integer.

    OUTPUT:

        The coefficient a_q(l1, l2, m1, m2, l) of Proposition 3.1. [q-KvPR].
    """
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


def calculate_recurrence(max_N, filename, num_cores=1):
    r"""
    calculate_recurrence calculates the recurrence relation for l=0 for the first 
    polynomials.

    INPUT:

        ``max_N``     - Integer, where max_N is the number of recurrence 
                        relations this function will compute.
        ``filename``  - String, filename where the recurrence relation is saved.
        ``num_cores`` - Positive Integer, by default 1, gives the number of 
                        cores to use for the computation.

    OUTPUT:

        A list ``recurrence`` with tuples as elements of the form 
        ``(a_{l-1/2, l-1/2}, a_{l+1/2, l+1/2}``.
    """

    print "Calculate the %i few recurrence relation elements (could take a while):" % MAX_N
    
    @parallel(num_cores) 
    def recurrence_a_parallel(l):
        if l == 0:
            return (0, recrel_a(q, 0, 0, 1/2, 1/2, 0))
        a1 = recrel_a(q, l/2, l/2, l/2-1/2, l/2-1/2, 0)
        a2 = recrel_a(q, l/2, l/2, l/2+1/2, l/2+1/2, 0)
        return (a1, a2)
    
    gen = recurrence_a_parallel(range(max_N))
    data = {}
    for rec in gen:
        print rec[0][0][0]
        data[rec[0][0][0]] = rec[1]
        save(data, REC_FILE)
    q = var('q')

def generate_weird_polynomials(recurrence):
    polynomials = [1, x/recurrence[0][1]**2]

    for i in range(2, len(recurrence)):
        newpoly = x/recurrence[i][1]**2*polynomials[i-1] \
            - (recurrence[i][0]/recurrence[i][1])**2*polynomials[i-2]
        polynomials.append(newpoly)
        print i,

    return polynomials

# Functional programming alert!
datafiles = ['data/rec10.sobj', 'data/rec20.sobj', 'data/rec30.sobj']
recurrence_relation = sum(map(lambda s : load(s), datafiles), [])

