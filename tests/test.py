import ctypes
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('dll')
parser.add_argument('mu0', type=float)
parser.add_argument('LWP', type=float)
parser.add_argument('r_e', type=float)
args = parser.parse_args()

testlib = ctypes.CDLL(args.dll)
run = testlib.run_spectral_test
run.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float]

run(args.mu0, args.LWP, args.r_e)
