#!/usr/bin/env python

import sys
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from math import log

# Read command line argument

scenario = sys.argv[1]

l2 = []
mesh = ['20', '40', '80', '160', '320']
h = [20, 40, 80, 160, 320]

# Read FD outputs
for i in range(0, len(mesh)):
	# Read second order values
	fname = './solns/soln'+scenario+mesh[i]+'.dat'
	f = open(fname, 'r')
	lines = f.read()

	match = re.findall('L2 Error Norm:.*', lines)
	match = match[0].split()
	val = match[-1]
	val = float(val)
	l2.append(val)
	f.close()
	
# Compute slope
err_0 = l2[0]
err_0 = log(err_0)
err_1 = l2[-1]
err_1 = log(err_1)
n_0 = log(h[0])
n_1 = log(h[-1])

slope = (err_1 - err_0)/(n_1 - n_0)

if scenario == '1d2nd':
	titlestr = 'Convergence Analysis Plot for 1D 2nd Order Problem'
elif scenario == '1d4th':
	titlestr = 'Convergence Analysis Plot for 1D 4th Order Problem'
elif scenario == '2d2nd':
	titlestr = 'Convergence Analysis Plot for 2D 2nd Order Problem'
elif scenarion == '2nd4th':
	titlestr = 'Convergence Analysis Plot for 2D 4th Order Problem'
else:
	raise ValueError('Invalid input')

# Plot
plt.loglog(h, l2, 'o-', color='g', label='slope={}'.format(round(slope, 3)))
plt.xlabel('# of Nodes')
plt.ylabel('L2 Error Norm')
plt.title(titlestr)
plt.legend(loc='upper right')
plt.savefig('./docs/'+scenario+'.png')

# Echo results to stdout
for i in range(0, len(h)):
	print('h: {} l2: {}'.format(h[i], l2[i]))

