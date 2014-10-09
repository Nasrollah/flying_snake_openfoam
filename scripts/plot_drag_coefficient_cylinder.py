#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_drag_coefficient_cylinder.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot drag coefficient over a cylinder and compare
#			   with Koumoutsakos and Leonard (1995).


import os
import argparse

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plot drag coefficient over a '
									 'cylinder and compare with results from '
									 'Koumoutsakos and Leonard (1995)',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill parser with arguments
	parser.add_argument('--case', dest='case_directory', type=str, default='.',
						help='directory of the simulation case')
	parser.add_argument('--start', dest='start', type=float, default=0.0,
						help='starting-time')
	parser.add_argument('--end', dest='end', type=float, default=3.5,
						help='ending-time')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the figure')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure')
	parser.add_argument('--compare', dest='compare', type=str,
						default='/home/mesnardo/flying_snake_openfoam/resources'
								'/cylinder_drag_coefficient_Re550_'
								'koumoutsakos_leonard_1995.dat',
						help='path of results from '
							 'Koumoutsakos and Leonard (1995)')
	parser.add_argument('--output', dest='output_name', type=str, 
						default='drag_coefficients',
						help='name of the output .png file')
	parser.set_defaults(save=True)
	return parser.parse_args()


def main():
	"""Plots the instantaneous drag coefficient
	and compares with Koumoutsakos and Leonard (1995).
	"""
	# parse the command-line
	args = read_inputs()

	# initialization
	t, cd = numpy.empty(0), numpy.empty(0)

	# read instantaneous drag coefficient from files
	force_path = '%s/postProcessing/forceCoeffs/' % args.case_directory
	folders = sorted(os.listdir(force_path))
	for folder in folders:
		with open(force_path+folder+'/forceCoeffs.dat') as infile:
			t_tmp, cd_tmp = numpy.loadtxt(infile, dtype=float, 
										  usecols=(0,2), delimiter='\t', 
										  unpack=True)
		t = numpy.append(t, t_tmp)
		cd = numpy.append(cd, cd_tmp)
	
	# calculate ending indices
	t_start = (t[0] if not args.start else args.start)
	t_end = (t[-1] if not args.end else args.end)
	i_start = numpy.where(t >= t_start)[0][0]
	i_end = numpy.where (t >= t_end)[0][0]-1

	# keep useful slices
	t = t[i_start:i_end].copy()
	cd = cd[i_start:i_end].copy()

	# read drag coefficient from Koumoutsakos and Leonard (1995)
	with open('%s' % args.compare, 'r') as infile:
		t_compare, cd_compare = numpy.loadtxt(infile, dtype=float, 
											  delimiter='\t', unpack=True)
		# conversion to same non-dimensional time
		t_compare *= 0.5

	# plot instantaneous drag coefficient
	pyplot.figure(figsize=(8, 8))
	pyplot.grid(True)
	pyplot.xlabel('Time', fontsize=16)
	pyplot.ylabel(r'$C_d$', fontsize=18)
	pyplot.plot(t, cd, label='OpenFOAM', color='b', ls='-', lw=2)
	pyplot.plot(t_compare, cd_compare, label='Koumoutsakos and Leonard (1995)', 
				color='r', lw=0, marker='o', markersize=6)
	pyplot.ylim(0.5, 2.0)
	pyplot.legend(loc='best', prop={'size': 14})

	# save the figure
	if args.save:
		# create images folder if not existing
		images_directory = '%s/images' % args.case_directory
		if not os.path.isdir(images_directory):
			os.makedirs(images_directory)
		pyplot.savefig('%s/%s.png' % (images_directory, args.output_name))

	# display the figure
	if args.show:
		pyplot.show()


if __name__ == '__main__':
	main()
