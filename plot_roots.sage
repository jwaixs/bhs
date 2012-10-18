roots = load('data/roots')

q_values = [1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8]
max_N = 20

for newq in q_values:
    points = []
    for i in range(1, max_N):
        zeros = roots[newq][i]
        points += [point(zip(zeros, len(zeros)*[i/max_N]), rgbcolor=hue(i/max_N))]
    sum(points).save('plots/roots_%s.png' % RR(newq))
