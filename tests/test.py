import ctypes
import argparse
import numpy
import numpy.ctypeslib

parser = argparse.ArgumentParser()
parser.add_argument('dll')
parser.add_argument('mu0', type=float)
parser.add_argument('LWP', type=float)
parser.add_argument('r_e', type=float)
args = parser.parse_args()

testlib = ctypes.CDLL(args.dll)
run = testlib.run_spectral_test
run.argtypes = [ctypes.c_double]*3 + [numpy.ctypeslib.ndpointer(dtype=numpy.double, shape=(24,), flags='CONTIGUOUS')]*5
T_DB = numpy.empty((24,))
T_DIF = numpy.empty((24,))
T_DIR = numpy.empty((24,))
R_DIF = numpy.empty((24,))
R_DIR = numpy.empty((24,))

run(args.mu0, args.LWP, args.r_e, T_DB, T_DIF, T_DIR, R_DIF, R_DIR)

from matplotlib import pyplot
fig = pyplot.figure()
ax = fig.add_subplot(2,1,1)
ax.plot(T_DB, label='T_DB: direct beam')
ax.plot(T_DIF, label='T_DIF: diffuse transmissivity for diffuse incident')
ax.plot(T_DIR, label='T_DIR: diffuse transmissivity for direct incident')
ax.legend()
ax.set_title('transmissivities')
ax = fig.add_subplot(2,1,2)
ax.plot(R_DIF, label='R_DIF: diffuse reflectivity for diffuse incident')
ax.plot(R_DIR, label='R_DIR: diffuse reflectivity for direct incident')
ax.legend()
ax.set_title('reflectivities')
pyplot.show()