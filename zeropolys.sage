q = var('q')
recurrence = load('data/rec10') + load('data/rec20')
polynomials = [1, x/recurrence[0][1]]

for i in range(1, 20):
    newpoly = x/recurrence[i][1]*polynomials[i-1] \
        - recurrence[i][0]/recurrence[i][1]*polynomials[i-2]
    polynomials.append(newpoly)
    print i
