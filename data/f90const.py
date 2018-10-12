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
   return tuple([numpy.array(r) for r in result])

def write(path, **kwargs):
   with open(path, 'w') as f:
       for name in sorted(kwargs.keys(), key=str.lower):
           data = kwargs[name]
           if isinstance(data, int):
              f.write('integer, parameter :: %s = %i\n' % (name, data))
           else:
              f.write('real(rk), parameter :: %s(%i) = (/ &\n' % (name, len(data)))
              for i in range(0, len(data), 10):
                  f.write('      %s%s &\n' % (', '.join(['%.6e_rk' % v for v in data[i:i+10]]), ',' if i + 10 < len(data) else ''))
              f.write('   /)\n')

