lambda_mins, lambda_maxs, a_s, bs, es, fs, cs, ds = [], [], [], [], [], [], [], []
with open('slingo.dat', 'rU') as f:
    for i in range(3):
        f.readline()
    for l in f:
        items = l.rstrip('\n').split()
        lambda_min, lambda_max, a, b, e, f, c, d = map(float, items)
        lambda_mins.append(lambda_min)
        lambda_maxs.append(lambda_max)
        a_s.append(a)
        bs.append(b)
        cs.append(c)
        ds.append(d)
        es.append(e)
        fs.append(f)
n = len(lambda_mins)

def write(f, arr, name, public=False):
    f.write('real(rk), parameter%s :: %s(%i) = (/%s/)\n' % (', public' if public else '', name, len(arr), ', '.join(['%g_rk' % v for v in arr])))

with open('../../src/slingo_const.inc', 'w') as f:
    f.write('real(rk), parameter :: nslingo = %i\n' % len(lambda_mins))
    write(f, lambda_mins, 'lambda_min', public=True)
    write(f, lambda_maxs, 'lambda_max', public=True)
    write(f, a_s, 'slingo_a')
    write(f, bs, 'slingo_b')
    write(f, cs, 'slingo_c')
    write(f, ds, 'slingo_d')
    write(f, es, 'slingo_e')
    write(f, fs, 'slingo_f')
