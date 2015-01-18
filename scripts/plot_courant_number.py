#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_courant_number.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot the instantaneous maximum Courant number


import argparse
import os

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plot the residuals from '
												 'a simpleFoam simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--case', dest='case_directory', type=str,
						default=os.getcwd(),
						help='directory of the OpenFoam simulation')
	parser.add_argument('--times', dest='times', type=float, nargs='+',
					    default=[None, None],
						help='time-interval to plot the Courant number')
	parser.add_argument('--limits', dest='limits', type=float, nargs='+',
						default=[None, None],
						help='y-limits to plot')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the figure')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure')
	parser.set_defaults(save=True)
	return parser.parse_args()


def get_courant_number(courant_path):
	"""Reads and stores the instantaneous maximum Courant number.

	Arguments
	---------
	courant_path -- path of the file containing the maximum Courant number.
	
	Returns
	-------
	times -- array with the discrete time values.
	courant -- array with the instantaneous maximum Courant number.
	"""
	with open(courant_path, 'r') as infile:
		times, courant = numpy.loadtxt(infile, delimiter='\t', unpack=True)
	return times, courant


def main():
	"""Plots the instantaneous maximum Courant number."""
	# parse the command-line
	args = read_inputs()

	# store the maximum Courant number
	logs_directory = '%s/logs' % args.case_directory
	times, courant = get_courant_number('%s/CourantMax_0' % logs_directory)

	# create the figure
	pyplot.figure(figsize=(10, 6))
	pyplot.grid(True)
	pyplot.xlabel('times', fontsize=16)
	pyplot.ylabel('maximum Courant number', fontsize=16)
	pyplot.plot(times, courant, color='b', linestyle='-', linewidth=2)
	if args.times:
		pyplot.xlim(args.times[0], args.times[1])
	if args.limits:
		pyplot.ylim(args.limits[0], args.limits[1])

	# save the figure
	if args.save:
		# create images folder if not existing
		images_directory = '%s/images' % args.case_directory
		if not os.path.isdir(images_directory):
			os.makedirs(images_directory)
		pyplot.savefig('%s/maximum_courant_number.png' % images_directory)

	# display the figure
	if args.show:
		pyplot.show()


if __name__ == '__main__':
	main()
