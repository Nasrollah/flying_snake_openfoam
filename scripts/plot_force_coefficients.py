#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_forceCoeffs.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot aerodynamic coefficients of the flying snake


import os
import argparse

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the lift and drag '
									'coefficients for a given simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--case', dest='case_directory', type=str, default='.',
						help='directory of the OpenFOAM case')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the aerodynamic coefficients')
	parser.add_argument('--start', dest='start', type=float, default=None,
						help='starting-time to compute mean coefficients '
							 ' and to plot instantaneous values')
	parser.add_argument('--end', dest='end', type=float, default=None,
						help='ending-time to compute mean coefficients '
							 ' and to plot instantaneous values')
	parser.add_argument('--cuibm', dest='cuibm_path', type=str, default=None,
						help='path of cuIBM force coefficients for comparison')
	parser.add_argument('--output', '-o', dest='output_name', type=str, 
						default='force_coefficients',
						help='name of the file generated (no extension)')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure as a .png file')
	parser.set_defaults(save=True)
	return parser.parse_args()


def main():
	"""Plots aerodynamic coefficients."""
	# parse the command-line
	args = read_inputs()

	# initialization
	t, cd, cl = numpy.empty(0), numpy.empty(0), numpy.empty(0)
	
	# read force coefficients from files
	force_path = '%s/postProcessing/forceCoeffs/' % args.case_directory
	folders = sorted(os.listdir(force_path))
	for folder in folders:
		with open(force_path+folder+'/forceCoeffs.dat') as infile:
			t_tmp, cd_tmp, cl_tmp = numpy.loadtxt(infile, dtype=float, 
											   	  usecols=(0,2,3), 
												  delimiter='\t', unpack=True)
		t = numpy.append(t, t_tmp)
		cd = numpy.append(cd, cd_tmp)
		cl = numpy.append(cl, cl_tmp)

	# calculates boundary indices
	t_start = (t[0] if not args.start else args.start)
	t_end = (t[-1] if not args.end else args.end)
	i_start = numpy.where(t >= t_start)[0][0]
	i_end = numpy.where(t >= t_end)[0][0]-1

	# keep useful slices
	t = t[i_start:i_end].copy()
	cd = cd[i_start:i_end].copy()
	cl = cl[i_start:i_end].copy()

	# compute mean coefficients
	cd_mean = cd.sum()/cd.size
	cl_mean = cl.sum()/cl.size

	# output print
	print '\t--> mean drag coefficient: %g' % cd_mean
	print '\t--> mean lift coefficient: %g' % cl_mean

	# plot instantaneous force coefficients
	pyplot.figure()
	pyplot.grid(True)
	pyplot.xlabel('Times', fontsize=16)
	pyplot.ylabel('Force Coefficients', fontsize=16)
	pyplot.plot(t, cd, label=r'$C_d$', color='r', ls='-', lw=2)
	pyplot.plot(t, cl, label=r'$C_l$', color='b', ls='-', lw=2)
	if args.cuibm_path:
		# read cuIBM force coefficients
		with open(args.cuibm_path, 'r') as infile:
			t_cuibm, cd_cuibm, cl_cuibm = numpy.loadtxt(infile, dtype=float, 
													 	delimiter='\t', 
													 	unpack=True)
		# keep useful slices
		i_start = numpy.where(t_cuibm >= t_start)[0][0]
		if t_end > t_cuibm[-1]:
			i_end = t_cuibm.size-1
		else:
			i_end = numpy.where(t_cuibm >= t_end)[0][0]-1
		t_cuibm = t_cuibm[i_start:i_end].copy()
		# multiplied by factor 2.0 because forgot in cuIBM
		cd_cuibm = 2.*cd_cuibm[i_start:i_end].copy()
		cl_cuibm = 2.*cl_cuibm[i_start:i_end].copy()
		
		# compute mean coefficients
		cd_cuibm_mean = cd_cuibm.sum()/cd_cuibm.size
		cl_cuibm_mean = cl_cuibm.sum()/cl_cuibm.size

		# output print
		print '\t--> mean drag coefficient with cuIBM: %g' % cd_cuibm_mean
		print '\t--> mean lift coefficient with cuIBM: %g' % cl_cuibm_mean

		# add cuIBM instantaneous values to the figure
		pyplot.plot(t_cuibm, cd_cuibm, label=r'$C_d$ - cuIBM', 
					color='r', ls='--')
		pyplot.plot(t_cuibm, cl_cuibm, label=r'$C_l$ - cuIBM', 
					color='b', ls='--')

	pyplot.legend(loc='best', prop={'size': 16})
	
	if args.save:
		# create images folder is not existing
		images_path = '%s/images' % args.case_directory
		if not os.path.isdir(images_path):
			os.makedirs(images_path)
		pyplot.savefig('%s/%s.png' % (images_path, args.output_name))
	
	if args.show:
		pyplot.show()


if __name__ == '__main__':
	main()
