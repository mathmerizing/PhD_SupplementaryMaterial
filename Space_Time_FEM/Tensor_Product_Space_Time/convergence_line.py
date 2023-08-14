import os
import argparse
import numpy as np

# dimension of FEM computations
DIM = 2

# parse flags to python script
# arguments are N1, error1, N2, error2, shift
parser = argparse.ArgumentParser(description='Create convergence line from 2 data points')
parser.add_argument('N1', metavar='N1', type=int, help='N1')
parser.add_argument('error1', metavar='error1', type=float, help='error1')
parser.add_argument('N2', metavar='N2', type=int, help='N2')
parser.add_argument('error2', metavar='error2', type=float, help='error2')
parser.add_argument('shift', metavar='shift', type=float, help='shift')
args = parser.parse_args()

# create data points
N1 = args.N1
error1 = args.error1
N2 = args.N2
error2 = args.error2

# get convergence rate as log(error)/log(N)
_rate = (np.log(error2) - np.log(error1))/(np.log(N2) - np.log(N1))
rounded_rate = round(-DIM * _rate)
print(f"Rate:  {-DIM * _rate} (~{rounded_rate})")

# compute a new point on the line going through (N1, error1+shift) with slope rounded_rate at N2
error3 = np.exp(np.log(error1+args.shift) - (rounded_rate / DIM) * (np.log(N2) - np.log(N1)))
print(f"({N1},{error1+args.shift})({N2},{error3})")