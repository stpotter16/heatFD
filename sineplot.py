#!/usr/bin/env python

import sys
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from math import log
import numpy as np

# Read command line argument

fname = sys.argv[1]

f = open(fname, 'r')

lines = f.readlines()

x = []
T = []

for i in range(25, len(lines) - 2):
	line = lines[i]
	split = line.split()
	x.append(float(split[1]))
	T.append(float(split[2]))

f.close()

Ax = 3 * np.pi
xx = np.linspace(0, 1, 100)
yy = np.cos(Ax * xx)

plt.plot(x, T, 'o-', color='g', label='Finite difference solution')
plt.plot(xx, yy, '-', color='r', label='Analytical solution')
plt.xlabel('x')
plt.ylabel('T')
plt.title('Finite difference vs Analytical Solution')
plt.legend(loc='upper right')
plt.savefig('./solns/solnplt.png')

