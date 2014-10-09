#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_force_coefficients.py
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
	parser.add_argument('--case', dest='case_directory', type=str, 
						default=os.getcwd(),
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
	parser.add_argument('--name', dest='image_name', type=str, 
						default='force_coefficients',
						help='name of the file generated (no extension)')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure as a .png file')
	parser.add_argument('--compare', dest='other_cases', type=str, nargs='+',
						help='directories of other cases for comparison')
	parser.set_defaults(save=True)
	return parser.parse_args()


class Case(object):
	"""Contains the path of the simulation."""
	def __init__(self, path):
		"""Stores the path of the simulation.
		
		Arguments
		---------
		path -- path of the simulation.
		"""
		self.path = path

	def get_mean_coefficients(self):
		"""Computes the mean force coefficients."""
		self.cd_mean = self.cd.sum()/self.cd.size
		self.cl_mean = self.cl.sum()/self.cl.size


class OpenFoamCase(Case):
	"""Contains force coefficients of an OpenFoam simulation."""
	def __init__(self, directory, t_start=None, t_end=None):
		"""Reads the force coefficients and computes the mean coefficients 
		between two ending-times.
		
		Arguments
		---------
		directory -- directory of the simulation.
		t_start, t_end -- boundary times (default: None, None).
		"""
		Case.__init__(self, directory)
		# read force coefficients
		self.read_coefficients(t_start, t_end)
		# compute mean coefficients
		self.get_mean_coefficients()

	def read_coefficients(self, t_start, t_end):
		"""Reads force coefficients from files.
		
		Arguments
		---------
		t_start, t_end -- boundary times.
		"""
		# initialization
		self.t = numpy.empty(0)
		self.cd, self.cl = numpy.empty(0), numpy.empty(0)
		# read force coefficients files and append to arrays
		forces_directory = '%s/postProcessing/forceCoeffs' % self.path
		folders = sorted(os.listdir(forces_directory))
		for folder in folders:
			force_path = '%s/%s/forceCoeffs.dat' % (forces_directory, folder)
			with open(force_path, 'r') as infile:
				t_tmp, cd_tmp, cl_tmp = numpy.loadtxt(infile, dtype=float, 
													  usecols=(0, 2, 3), 
													  delimiter='\t', 
													  unpack=True)
			self.t = numpy.append(self.t, t_tmp)
			self.cd = numpy.append(self.cd, cd_tmp)
			self.cl = numpy.append(self.cl, cl_tmp)
		# calculate boundary indices
		t_start = (self.t[0] if not t_start else t_start)
		t_end = (self.t[-1] if not t_end else t_end)
		i_start = numpy.where(self.t >= t_start)[0][0]
		i_end = numpy.where(self.t >= t_end)[0][0]-1
		# keep useful slices
		self.t = self.t[i_start:i_end].copy()
		self.cd = self.cd[i_start:i_end].copy()
		self.cl = self.cl[i_start:i_end].copy()


class CuIBMCase(Case):
	"""Contains force coefficients of a cuIBM  simulation."""
	def __init__(self, path, t_start=None, t_end=None):
		"""Reads the force coefficients and computes the mean coefficients 
		between two ending-times.
		
		Arguments
		---------
		path -- path of the cuIBM force coefficients file.
		t_start, t_end -- boundary times (default: None, None).
		"""
		Case.__init__(self, path)
		# read force coefficients
		self.read_coefficients(t_start, t_end)
		# compute mean coefficients
		self.get_mean_coefficients()

	def read_coefficients(self, t_start, t_end):
		"""Reads force coefficients from cuIBM results.
		
		Arguments
		---------
		t_start, t_end -- boundary times.
		"""
		# read the file
		with open(self.path, 'r') as infile:
			self.t, self.cd, self.cl = numpy.loadtxt(infile, dtype=float, 
													 delimiter='\t', 
													 unpack=True)
		# keep useful slices
		i_start = numpy.where(self.t >= t_start)[0][0]
		if t_end > self.t[-1]:
			i_end = self.t.size-1
		else:
			i_end = numpy.where(self.t >= t_end)[0][0]-1
		self.t = self.t[i_start:i_end].copy()
		# multiply by factor 2.0 which does not appear in cuIBM
		self.cd = 2.*self.cd[i_start:i_end].copy()
		self.cl = 2.*self.cl[i_start:i_end].copy()
		

def plot_coefficients(cases, image_name=None, save=False, show=False):
	"""Plots force coefficients from different simulations.
	
	Arguments
	--------
	cases -- dictionary that contains info about all simulations to plot.
	image_name -- name of the file generated (default: None).
	save -- does save the figure as a .png file (default: False).
	show -- dislays the figure (default: False).
	"""
	# figure parameters
	pyplot.figure(figsize=(10, 6))
	pyplot.grid(True)
	pyplot.xlabel('time', fontsize=16)
	pyplot.ylabel('force coefficients', fontsize=16)
	# plot the main OpenFoam force coefficients
	pyplot.plot(cases['main'].t, cases['main'].cd, 
				label=r'$C_d$ - %s' % os.path.basename(cases['main'].path),
				color='r', ls='-', lw=2)
	pyplot.plot(cases['main'].t, cases['main'].cl, 
				label=r'$C_l$ - %s' % os.path.basename(cases['main'].path),
				color='b', ls='-', lw=2)
	# plot other OpenFoam force coefficients
	colors = ['g', 'c', 'm', 'y']
	for i, case in enumerate(cases['others']):
		pyplot.plot(case.t, case.cd,
					label=r'$C_d$ - %s' % os.path.basename(case.path),
					color=colors[i], ls='-', lw=1)
		pyplot.plot(case.t, case.cl,
					label=r'$C_l$ - %s' % os.path.basename(case.path),
					color=colors[i], ls='--', lw=1)
	# plot cuIBM force coefficients
	if cases['cuibm']:
		pyplot.plot(cases['cuibm'].t, cases['cuibm'].cd,
					label=r'$C_d$ - cuIBM', color='k', ls='-', lw=1)
		pyplot.plot(cases['cuibm'].t, cases['cuibm'].cl,
					label=r'$C_l$ - cuIBM', color='k', ls='--', lw=1)
	pyplot.legend(loc='upper left', prop={'size': 'small'},
						bbox_to_anchor=(1.0,1.0))
	# save the figure as a .png file
	if save:
		# create images folder if not existing
		images_directory = '%s/images' % cases['main'].path
		if not os.path.isdir(images_directory):
			os.makedirs(images_directory)
		pyplot.savefig('%s/%s.png' % (images_directory, image_name),
					   bbox_inches='tight')
	# display the figure
	if show:
		pyplot.show()


def main():
	"""Plots aerodynamic coefficients."""
	# parse the command-line
	args = read_inputs()

	# store case directories
	cases = {'main': args.case_directory,
			 'others': (None if not args.other_cases else args.other_cases),
			 'cuibm': (None if not args.cuibm_path else args.cuibm_path)}

	# read coefficients and compute mean values of main OpenFoam simulation
	cases['main'] = OpenFoamCase(cases['main'], 
								 t_start=args.start, t_end=args.end)

	# read coefficients and compute mean values of other OpenFoam simulations
	for i, directory in enumerate(cases['others']):
		cases['others'][i] = OpenFoamCase(directory, 
										  t_start=args.start, t_end=args.end)

	# read coefficients and compute mean values of cuIBM simulation
	if args.cuibm_path:
		cases['cuibm'] = CuIBMCase(cases['cuibm'], 
								   t_start=args.start, t_end=args.end)

	# plot force coefficients
	plot_coefficients(cases, 
					  image_name=args.image_name, 
					  save=args.save, 
					  show=args.show)


if __name__ == '__main__':
	main()
