def cf2elm(cf):
    if len(cf) == 0:
        return 0
    elif len(cf) == 1:
        return cf[0]
    else:
        tail = cf2elm(cf[1:])
        if tail == 0:
            return cf[0]
        else:
            return cf[0] + 1/tail

from itertools import product

elements = [-4, -2*I-3, -I-3, -3, I-3, 2*I-3, -3*I-2, -2*I-2, -I-2, -2, \
    I-2, 2*I-2, 3*I-2, -3*I-1, -2*I-1, -I-1, -1, I-1, 2*I-1, 3*I-1, -4*I, \
    -3*I, -2*I, -I, I, 2*I, 3*I, 4*I, -3*I+1, -2*I+1, -I+1, 1, I+1, 2*I+1, \
    3*I+1, -3*I+2, -2*I+2, -I+2, 2, I+2, 2*I+2, 3*I+2, -2*I+3, -I+3, 3, \
    I+3, 2*I+3, 4]

ccf4_approx = product(elements, repeat=1)

points = []
for (cfa, cfb, cfc) in product(ccf4_approx, repeat=3):
    print cfa, cfb, cfc
    cor = cf2elm(cfa) + cf2elm(cfb) + cf2elm(cfc)
    points.append(point((real_part(cor), imag_part(cor)), size=1))
sum(points).save('test.png')
