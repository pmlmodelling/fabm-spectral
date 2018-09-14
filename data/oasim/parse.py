import numpy

def parse(path, nskip, columns):
   result = [[] for i in range(len(columns))]
   with open(path, 'rU') as f:
       for i in range(nskip):
           f.readline()
       for l in f:
           items = l.rstrip('\n').split()
           for i, r in zip(columns, result):
              r.append(float(items[i]))
   return tuple(result)

def write(path, **kwargs):
   with open(path, 'w') as f:
       for name, data in kwargs.items():
           if isinstance(data, int):
              f.write('integer, parameter :: %s = %i\n' % (name, data))
           else:
              f.write('real(rk), parameter :: %s(%i) = (/%s/)\n' % (name, len(data), ', '.join(['%e_rk' % v for v in data])))

lam, a, b = parse('abw25b.dat', 6, range(3))
write('../../src/water_const.inc', nlambda_w=len(lam), lambda_w=lam, a_w=a, b_w=b)

lambda_min, lambda_max, a, b, e, f, c, d = parse('slingo.dat', 3, range(8))
write('../../src/slingo_const.inc', nlambda_slingo=len(lambda_min), lambda_min_slingo=lambda_min, lambda_max_slingo=lambda_max, a_slingo=a, b_slingo=b, c_slingo=c, d_slingo=d, e_slingo=e, f_slingo=f)
