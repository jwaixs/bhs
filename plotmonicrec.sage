from sage.symbolic.expression_conversions import RingConverter

q = var('q')
data = load('data/rec10') #+ load('data/rec20') + load('data/rec30') \
    #+ load('data/rec40') + load('data/rec50')
CORES = 7

def plot_function_list(f, begin=0.01, end=1.0, steps=0.005):
    result = []
    for r in srange(begin, end, steps):
        result.append((r, f.substitute(q=r)))
    return result

def plot_monic_rec(n, color='blue'):
    coef = data[n][0]**2 * data[n-1][1]**2
    return list_plot(plot_function_list(coef), plotjoined=True, rgbcolor=color)

plots = []
for i in range(1, 10):
    print 'plotting %i' % i
    plots.append(plot_monic_rec(i, color=hue(i/10.0)))
