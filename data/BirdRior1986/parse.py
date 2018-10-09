lams, exters, a_ws, a_os, a_us = [], [], [], [], []
with open('table1.txt', 'rU') as f:
    labels = f.readline().rstrip('\n').split('\t')
    for l in f:
        items = l.rstrip('\n').split('\t')
        lam, exter, a_w, a_o, a_u = map(float, items)
        lams.append(lam*1000)
        exters.append(exter)
        a_ws.append(a_w)
        a_os.append(a_o)
        a_us.append(a_u)

import sys
sys.path.append('..')
import f90const
f90const.write('../../src/birdrior1986_const.inc', nlambda_birdrior1986=len(lams), lambda_birdrior1986=lams, a_w_birdrior1986=a_ws, a_o_birdrior1986=a_os, a_u_birdrior1986=a_us)

from matplotlib import pyplot
fig = pyplot.figure()
ax = fig.gca()
ax.plot(lams, a_us, '-')
ax.grid(True)
pyplot.show()