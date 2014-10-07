#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_cylinder_data_file.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate cylinder coordinates and write in a data file


import os
import argparse
import math

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates a data file '
									 'eith 2D-cylinder coordinates',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('-n', dest='n', type=int, default=100,
						help='number of points on the cylinder')
	parser.add_argument('--center', '-c', dest='center', type=float, nargs='+',
						default=[0.0, 0.0],
						help='coordinates of the center of the cylinder')
	parser.add_argument('--radius', '-r', dest='r', type=float, default=0.5,
						help='radius of the cylinder')
	parser.add_argument('--save-dir', '-s', dest='save_dir', type=str,
						default='.',
						help='directory where the file will be saved')
	parser.add_argument('--show', dest='show', action='store_true',
						help='plots and displays the cylinder in a figure')
	return parser.parse_args()


def main():
	"""Generates cylinder coordinates and writes into a data file."""
	# parse the command-line
	args = read_inputs()

	# store parameters
	R = args.r
	N = args.n
	x_center, y_center = args.center[0], args.center[1]

	# create the cylinder points
	x = x_center + R*numpy.cos(numpy.linspace(0.0, 2.0*math.pi, N+1)[:-1])
	y = y_center + R*numpy.sin(numpy.linspace(0.0, 2.0*math.pi, N+1)[:-1])

	# write coordinates into a data file
	with open('%s/cylinder_%g_%d.dat' % (args.save_dir, R, N), 'w') as outfile:
		outfile.write('%d\n' % N)
		numpy.savetxt(outfile, numpy.c_[x, y], fmt='%.6f', delimiter='\t')

	# plot and display the cylinder
	if args.show:
		pyplot.figure(figsize=(6, 6))
		pyplot.grid(True)
		pyplot.xlabel(r'$x$', fontsize=18)
		pyplot.ylabel(r'$y$', fontsize=18)
		pyplot.plot(numpy.append(x, x[0]), numpy.append(y, y[0]),
					color='b', ls='-', lw=2, marker='o', markersize=6)
		pyplot.title('cylinder (%d points)' % N)
		pyplot.show()


if __name__ == '__main__':
	main()
