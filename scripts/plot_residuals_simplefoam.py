#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_residuals_simplefoam.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot the residuals from a simpleFoam simulation


import argparse

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plot the residuals from '
												 'a simpleFoam simulation')
	# fill the parser with arguments
	parser.add_argument('--case', dest='case_directory', type=str,
						help='directory of the simpleFoam simulation')
	parser.add_argument('--show', dest='show', action='store_true',
						help='display the figure')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='do not save the figure')
	parser.set_defaults(save=True)
	return parser.parse_args()


def get_residuals(residuals_path):
	"""Reads and stores the residuals.

	Arguments
	---------
	residuals_path -- path of the file containing the reisudals.
	
	Returns
	-------
	ites -- array with the iteration numbers.
	residuals -- array with the residuals.
	"""
	with open(residuals_path, 'r') as infile:
		ites, residuals = numpy.loadtxt(infile, delimiter='\t', unpack=True)
	return ites, residuals


def main():
	"""Plots the residuals obtained in a simpleFoam simulation."""
	# parse the command-line
	args = read_inputs()

	# store the residuals
	residuals_directory = '%s/logs' % args.case_directory
	ites_u, residuals_u = get_residuals('%s/Ux_0' % residuals_direcotry)
	ites_v, residuals_v = get_residuals('%s/Uy_0' % residuals_direcotry)
	ites_p, residuals_p = get_residuals('%s/p_0' % residuals_direcotry)

	# create the figure
	pyplot.figure(figsize=(8, 8))
	pyplot.grid(True)
	pyplot.xlabel('iteration number', fontisze=16)
	pyplot.ylabel('reisudals', fontsize=16)
	pyplot.plot(ites_p, residuals_p, label=r'$p$', color='r', ls='-', lw=2)
	pyplot.plot(ites_u, residuals_u, label=r'$U_x$', color='b', ls='-', lw=2)
	pyplot.plot(ites_v, residuals_v, label=r'$U_y$', color='g', ls='-', lw=2)
	pyplot.legend(loc='best', prop={'size': 18})

	# save the figure
	if args.save:
		# create images folder if not existing
		images_directory = '%s/images' % args.case_directory
		if not os.path.isdir(images_directory):
			os.makedirs(images_directory)
		pyplot.savefig('%s/residuals.png' % images_directory)

	# display the figure
	if args.show:
		pyplot.show()


if __name__ == '__main__':
	main()
