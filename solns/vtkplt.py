#!/usr/bin/env python

import h5py

h = h5py.File('2dplot.h5', 'r')

nodedata = h.get('nodes')
nodearr = nodedata[:]
nodes = nodearr[0]

spacing = (1.0 - 0.0)/(nodes - 1)

tempdata = h.get('T')
temparr = tempdata[:]

total = nodes * nodes

f = open('2dplot.vtk', 'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('Temperature data from class project\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_POINTS\n')
f.write('DIMENSIONS {} {} 1 \n'.format(nodes, nodes))
f.write('SPACING {} {} 0.0000\n'.format(spacing, spacing))
f.write('ORIGIN 0.0000 0.0000 0.0000\n\n')

f.write('POINT_DATA {}\n'.format(total))
f.write('SCALARS temperature double 1\n')
f.write('LOOKUP_TABLE default\n')
for i in range(len(temparr)):
	f.write('{}\n'.format(temparr[i][0]))
	
f.close()