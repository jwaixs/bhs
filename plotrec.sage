q = var('q')
recurrence = load('recurrence')

recplots = []
for n in range(10):
    recplots.append(plot(recurrence[n][0], q, 0, 1, rgbcolor=hue(n/10.0)))
    recplots.append(plot(recurrence[n][1], q, 0, 1, rgbcolor=hue(n/10.0)))
    recplots.append(text(n/2, (0.1, 1-n/20.0), rgbcolor=hue(n/10.0)))
    recplots.append(text(n/2, (0.1, n/20.0), rgbcolor=hue(n/10.0)))

for n in range(10):
    newplot = []
    newplot.append(plot(recurrence[n][0], q, 0, 1, rgbcolor=hue(n/10.0)))
    newplot.append(plot(recurrence[n][1], q, 0, 1, rgbcolor=hue(n/10.0)))
    newplot.append(text('0:%s' % str(n/2), (0.1, recurrence[n][0].substitute(q=0.1))))
    newplot.append(text('1:%s' % str(n/2), (0.1, recurrence[n][1].substitute(q=0.1))))
    sum(newplot).save('recurrenceplot_%i.png' % n)

newplot = []
for n in range(14):
    newplot.append(point((n, recurrence[n][0].substitute(q=0.5)), color="pink"))
    newplot.append(point((n, recurrence[n][1].substitute(q=0.5))))

sum(newplot).save('recurrencepoints.png')
