#!/usr/bin/env python

# file: $OPENFOAM_CASE_DIR/postProcessing/plot_forceCoeffs.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: plots aerodynamic coefficients of the flying snake


import os
import argparse

import numpy as np
from matplotlib import pyplot as plt


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the lift and drag '
									'coefficients for a given simulation')
	# fill the parser with arguments
	parser.add_argument('--case', dest='case', type=str, default='.',
						help='path of the OpenFOAM case')
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
	parser.add_argument('--output', '-o', dest='output', type=str, 
						default='force_coefficients',
						help='name of the file generated')
	return parser.parse_args()


def main():
	"""Plots aerodynamic coefficients."""
	# parse the command-line
	args = read_inputs()

	# initialization
	t, cd, cl = np.empty(0), np.empty(0), np.empty(0)
	
	# read force coefficients from files
	force_path = '%s/postProcessing/forceCoeffs/' % args.case
	folders = sorted(os.listdir(force_path))
	for folder in folders:
		with open(force_path+folder+'/forceCoeffs.dat') as infile:
			t_tmp, cd_tmp, cl_tmp = np.loadtxt(infile, dtype=float, 
											   usecols=(0,2,3), delimiter='\t',
											   unpack=True)
		t = np.append(t, t_tmp)
		cd = np.append(cd, cd_tmp)
		cl = np.append(cl, cl_tmp)

	# calculates boundary indices
	t_start = (t[0] if not args.start else args.start)
	t_end = (t[-1] if not args.end else args.end)
	i_start = np.where(t >= t_start)[0][0]
	i_end = np.where(t >= t_end)[0][0]-1

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
	plt.figure()
	plt.grid(True)
	plt.xlabel('Times', fontsize=16)
	plt.ylabel('Force Coefficients', fontsize=16)
	plt.plot(t, cd, label=r'$C_d$', color='r', ls='-', lw=2)
	plt.plot(t, cl, label=r'$C_l$', color='b', ls='-', lw=2)
	if args.cuibm_path:
		# read cuIBM force coefficients
		with open(args.cuibm_path, 'r') as infile:
			t_cuibm, cd_cuibm, cl_cuibm = np.loadtxt(infile, dtype=float, 
													 delimiter='\t', 
													 unpack=True)
		# keep useful slices
		i_start = np.where(t_cuibm >= t_start)[0][0]
		i_end = np.where(t_cuibm >= t_end)[0][0]-1
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
		plt.plot(t_cuibm, cd_cuibm, label=r'$C_d$ - cuIBM', color='r', ls='--')
		plt.plot(t_cuibm, cl_cuibm, label=r'$C_l$ - cuIBM', color='b', ls='--')

	plt.legend(loc='best', prop={'size': 16})
	plt.savefig('%s/postProcessing/%s.png' % (args.case, args.output))
	if args.show:
		plt.show()


if __name__ == '__main__':
	main()
