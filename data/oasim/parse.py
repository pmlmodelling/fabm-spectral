import numpy
import sys
sys.path.append('..')
import f90const

def parse(path, nskip, columns):
   result = [[] for i in range(len(columns))]
   with open(path, 'rU') as f:
       for i in range(nskip):
           f.readline()
       for l in f:
           items = l.rstrip('\n').split()
           for i, r in zip(columns, result):
              r.append(float(items[i]))
   return tuple([numpy.array(r) for r in result])

lam, a, b = parse('abw25b.dat', 6, range(3))
f90const.write('../../src/water_const.inc', nlambda_w=len(lam), lambda_w=lam, a_w=a, b_w=b)

lambda_min, lambda_max, a, b, e, f, c, d = parse('slingo.dat', 3, range(8))
a *= 0.01
f *= 0.001
lambda_min *= 1000
lambda_max *= 1000
f90const.write('../../src/slingo_const.inc', nlambda_slingo=len(lambda_min), lambda_min_slingo=lambda_min, lambda_max_slingo=lambda_max, lambda_slingo=(lambda_min + lambda_max)*0.5, a_slingo=a, b_slingo=b, c_slingo=c, d_slingo=d, e_slingo=e, f_slingo=f)

phyto_data = {}
with open('acbc25b.dat', 'rU') as f:
    current_data = None
    for l in f:
        items = l.rstrip('\n').split()
        if len(items) == 1:
            current_data = phyto_data.setdefault(items[0], [])
        elif current_data is not None:
            assert len(items) == 3
            current_data.append(map(float, items))
k2v = {}
for name, current_data in phyto_data.items():
    current_data = numpy.array(current_data)
    k2v['lambda_%s' % name] = current_data[:, 0]
    k2v['a_%s' % name] = current_data[:, 1]
    k2v['b_%s' % name] = current_data[:, 2]
f90const.write('../../src/phyto_const.inc', **k2v)

lam, ET, tau_r, a_o, a_v, a_o2, a_co2 = parse('atmo25b.dat', 4, range(7))
f90const.write('../../src/oasim_const.inc', nlambda_oasim=len(lam), lambda_oasim=lam, ET_oasim=ET, tau_r_oasim=tau_r, a_o_oasim=a_o, a_v_oasim=a_v, a_u_oasim=a_o2 + a_co2)

lams_br, a_us_br = [], []
with open('../birdrior1986/table1.txt', 'rU') as f:
    labels = f.readline().rstrip('\n').split('\t')
    for l in f:
        items = l.rstrip('\n').split('\t')
        lams_br.append(float(items[0])*1000)
        a_us_br.append(float(items[4]))

from matplotlib import pyplot
fig = pyplot.figure(figsize=(12,5))
ax = fig.add_subplot(121)
BR_tau_r = 1. / (115.6406 * (lam/1000)**4 - 1.335 * (lam/1000)**2)
ax.plot(lam, tau_r, '-', label='OASIM')
ax.plot(lam, BR_tau_r, '-', label='Bird & Riordan 1986')
ax.legend()
ax.grid(True)

ax = fig.add_subplot(122)
ax.plot(lam, a_o2 + a_co2, '-', label='OASIM')
ax.plot(lams_br, a_us_br, '-', label='Bird & Riordan 1986')
ax.legend()
ax.grid(True)

fig.savefig('comparison.png', dpi=300)
