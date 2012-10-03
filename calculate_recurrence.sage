load('bhs.sage')
load('recrel.sage')

NUM_CORES = 4
REC_FILE = 'recurrence'
MAX_N = 100

q = var('q')

@parallel(NUM_CORES) 
def recurrence_a_parallel(l):
    if l == 0:
        return (0, recrel_a(q, 0, 0, 1/2, 1/2, 0))
    a1 = recrel_a(q, l/2, l/2, l/2-1/2, l/2-1/2, 0)
    a2 = recrel_a(q, l/2, l/2, l/2+1/2, l/2+1/2, 0)
    return (a1, a2)


print "Calculate the %i few recurrence relation elements (could take a while):" % MAX_N
gen = recurrence_a_parallel(range(MAX_N))
data = {}
for rec in gen:
    print rec[0][0][0]
    data[rec[0][0][0]] = rec[1]
    save(data, REC_FILE)
