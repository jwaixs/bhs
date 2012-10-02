q = var('q')
recurrence = load('recurrence')

recplots = []
for n in range(10):
    recplots.append(plot(recurrence[n][0], q, 0, 1, rgbcolor=hue(n/10.0)))
    recplots.append(plot(recurrence[n][1], q, 0, 1, rgbcolor=hue(n/10.0)))
    recplots.append(text(n/2, (0.1, 1-n/20.0), rgbcolor=hue(n/10.0)))
    recplots.append(text(n/2, (0.1, n/20.0), rgbcolor=hue(n/10.0)))
sum(recplots)
