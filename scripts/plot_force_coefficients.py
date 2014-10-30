#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_force_coefficients.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot aerodynamic coefficients of an OpenFoam case


import os
import sys
import argparse
import logging
import datetime

import numpy
from scipy import signal
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
						help='directory of the OpenFoam case')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the aerodynamic coefficients')
	parser.add_argument('--times', dest='times', type=float, nargs='+', 
						default=[None, None],
						help='time-interval to compute mean coefficients')
	parser.add_argument('--limits', dest='limits', type=float, nargs='+',
						default=[None, None, None, None],
						help='limits of the plot')
	parser.add_argument('--legend', dest='legend', type=str, nargs='+',
						help='legend for each simulation to plot')
	parser.add_argument('--no-lift', dest='lift', action='store_false',
						help='does not plot the lift coefficients')
	parser.add_argument('--no-drag', dest='drag', action='store_false',
						help='does not plot the drag coefficients')
	parser.add_argument('--name', dest='image_name', type=str, 
						default='force_coefficients',
						help='name of the file generated (no extension)')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure as a .png file')
	parser.add_argument('--compare', dest='other_cases', type=str, nargs='+',
						help='directories of other cases for comparison')
	parser.add_argument('--cuibm', dest='cuibm_path', type=str, default=None,
						help='path of cuIBM force coefficients for comparison')
	parser.add_argument('--kl1995', dest='kl1995', action='store_true',
						help='plots instantaneous drag coefficient from '
							 'Koumoutsakos and Leonard (1995)')
	parser.set_defaults(save=True, lift=True, drag=True)
	# parse the command-line
	args = parser.parse_args()
	# log the command-line arguments
	if args.save:
		log_path = '%s/%s.log' % (args.case_directory, 
								  os.path.basename(__file__))
		logging.basicConfig(filename=log_path, format='%(message)s', 
							level=logging.INFO)
		logging.info('\n\n%s' % str(datetime.datetime.now()))
		logging.info(' '.join(sys.argv))
		logging.info(args)
	return args


class Case(object):
	"""Contains the path of the simulation."""
	def __init__(self, path):
		"""Stores the path of the simulation.
		
		Arguments
		---------
		path -- path of the simulation.
		"""
		self.path = path

	def get_time_limits(self, t_start, t_end):
		"""Computes the time-limits and their indices in the time array.
		
		Arguments
		---------
		t_start, t_end -- time-limits to compute mean coefficients and Strouhal.
		"""
		self.t_start = (self.t[0] if not t_start else t_start)
		self.t_end = (self.t[-1] if not t_end else t_end)
		self.i_start = numpy.where(self.t >= self.t_start)[0][0]
		self.i_end = numpy.where(self.t >= t_end)[0][0]-1

	def get_mean_coefficients(self):
		"""Computes the mean force coefficients."""
		self.cd_mean = (self.cd[self.i_start:self.i_end].sum()
						/ self.cd[self.i_start:self.i_end].size)
		self.cl_mean = (self.cl[self.i_start:self.i_end].sum()
						/ self.cl[self.i_start:self.i_end].size)

	def get_extremum_coefficients_old(self):
		"""Computes the extrema of the force coefficients."""
		self.cd_max = self.cd[self.i_start:self.i_end].max()
		self.cd_min = self.cd[self.i_start:self.i_end].min()
		self.cl_max = self.cl[self.i_start:self.i_end].max()
		self.cl_min = self.cl[self.i_start:self.i_end].min()

	def get_extremum_coefficients(self):
		"""Computes the extrema of the force coefficients 
		reached during the last period and computes the Strouhal number.
		"""
		def get_extremum_indices(x):
			"""Returns the index of the extrema within a given array.
			
			Arguments
			---------
			x -- array where extremum indices will be computed.

			Returns
			-------
			minima -- index of the minima.
			maxima -- index of the maxima.
			"""
			minima = signal.argrelextrema(x, numpy.less, order=5)[0]
			maxima = signal.argrelextrema(x, numpy.greater, order=5)[0]
			if minima.size == 0:
				minima = numpy.array([x.min()])
			if maxima.size == 0:
				maxima = numpy.array([x.max()])
			return minima, maxima
		
		# store useful slices
		t = self.t[self.i_start:self.i_end]
		cd = self.cd[self.i_start:self.i_end]
		cl = self.cl[self.i_start:self.i_end]
		# compute extremum drag coefficients
		minima, maxima = get_extremum_indices(cd)
		self.cd_min, self.cd_max = cd[minima[-1]], cd[maxima[-1]]
		# compute extremum lift coefficients
		minima, maxima = get_extremum_indices(cl)
		self.cl_min, self.cl_max = cl[minima[-1]], cl[maxima[-1]]
		# calculate the Strouhal number (chord-length=1.0 velocity=1.0)
		if minima.size > 1:
			self.strouhal = 1.0/(t[minima[-1]]-t[minima[-2]])
		else:
			self.strouhal = 0.0

	def get_strouhal_number_old(self):
		"""Calculates the Strouhal number."""
		spectrum = numpy.fft.fft(self.cl[self.i_start:self.i_end])
		dt = self.t[self.i_start+1] - self.t[self.i_start]
		frequencies = numpy.fft.fftfreq(spectrum.size, dt)
		mask = frequencies > 0.0
		peaks = frequencies[mask]
		magns = numpy.absolute(spectrum[mask])
		self.strouhal =  peaks[numpy.argmax(magns)]


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
		self.read_coefficients()
		# compute mean coefficients
		self.get_time_limits(t_start, t_end)
		self.get_mean_coefficients()
		# compute extrema and Strouhal number
		self.get_extremum_coefficients()

	def read_coefficients(self):
		"""Reads force coefficients from files."""
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
		self.read_coefficients()
		# compute mean coefficients
		self.get_time_limits(t_start, t_end)
		self.get_mean_coefficients()

	def read_coefficients(self):
		"""Reads force coefficients from cuIBM results."""
		# read the file
		with open(self.path, 'r') as infile:
			self.t, self.cd, self.cl = numpy.loadtxt(infile, dtype=float, 
													 delimiter='\t', 
													 unpack=True)
		# cuIBM script does not include the coefficient 2.0
		self.cd *= 2.0
		self.cl *= 2.0
		

def plot_coefficients(cases, args):
	"""Plots force coefficients from different simulations.
	
	Arguments
	--------
	cases -- dictionary that contains info about all simulations to plot.
	args -- namespace of command-line arguments
	"""
	# get the legend for each plot
	legend = {}
	if not args.legend:
		legend['main'] = os.path.basename(cases['main'].path)
		legend['others'] = []
		for case in cases['others']:
			legend['others'].append(os.path.basename(case.path))
	else:
		legend['main'] = args.legend[0]
		legend['others'] = []
		for i in xrange(len(cases['others'])):
			legend['others'].append(args.legend[i+1])
	# figure parameters
	pyplot.figure(figsize=(10, 6))
	pyplot.grid(True)
	pyplot.xlabel('time', fontsize=16)
	pyplot.ylabel('force coefficients', fontsize=16)
	# plot the main OpenFoam force coefficients
	if args.drag:
		pyplot.plot(cases['main'].t, cases['main'].cd, 
					label=r'$C_d$ - %s' % legend['main'], 
					color='r', ls='-', lw=2)
	if args.lift:
		pyplot.plot(cases['main'].t, cases['main'].cl, 
					label=r'$C_l$ - %s' % legend['main'], 
					color='b', ls='-', lw=2)
	# plot other OpenFoam force coefficients
	colors = ['g', 'c', 'm', 'y', 'k']
	for i, case in enumerate(cases['others']):
		if args.drag:
			pyplot.plot(case.t, case.cd, 
						label=r'$C_d$ - %s' % legend['others'][i],
						color=colors[i], ls='-', lw=1)
		if args.lift:
			pyplot.plot(case.t, case.cl,
						label=r'$C_l$ - %s' % legend['others'][i],
						color=colors[i], ls='--', lw=1)
	# plot cuIBM force coefficients
	if cases['cuibm']:
		if args.drag:
			pyplot.plot(cases['cuibm'].t, cases['cuibm'].cd,
						label=r'$C_d$ - cuIBM', color='k', ls='-', lw=1)
		if args.drag:
			pyplot.plot(cases['cuibm'].t, cases['cuibm'].cl,
						label=r'$C_l$ - cuIBM', color='k', ls='--', lw=1)
	# plot Koumoutsakos and Leonard (1995) drag coefficients
	if args.kl1995:
		pyplot.plot(cases['kl1995'].t, cases['kl1995'].cd,
					label=r'$C_d$ - Koumoutsakos and Leonard (1995)', 
					color='k', lw=0, marker='o', markersize=6)
	x_min = (cases['main'].t[0] if not args.limits else args.limits[0])
	x_max = (cases['main'].t[-1] if not args.limits else args.limits[1])
	y_min = (-2.0 if not args.limits else args.limits[2])
	y_max = (+2.0 if not args.limits else args.limits[3])
	pyplot.xlim(x_min, x_max)
	pyplot.ylim(y_min, y_max)
	pyplot.legend(loc='upper left', prop={'size': 'small'},
						bbox_to_anchor=(1.0,1.0))
	# save the figure as a .png file
	if args.save:
		# create images folder if not existing
		images_directory = '%s/images' % args.case_directory
		if not os.path.isdir(images_directory):
			os.makedirs(images_directory)
		pyplot.savefig('%s/%s.png' % (images_directory, args.image_name),
					   bbox_inches='tight')
	# display the figure
	if args.show:
		pyplot.show()


def main():
	"""Plots aerodynamic coefficients."""
	# parse the command-line
	args = read_inputs()

	# store case directories
	cases = {'main': args.case_directory,
			 'others': ([] if not args.other_cases else args.other_cases),
			 'cuibm': (None if not args.cuibm_path else args.cuibm_path)}

	# read coefficients and compute mean values of main OpenFoam simulation
	cases['main'] = OpenFoamCase(cases['main'], 
								 t_start=args.times[0], t_end=args.times[1])

	# read coefficients and compute mean values of other OpenFoam simulations
	for i, directory in enumerate(cases['others']):
		cases['others'][i] = OpenFoamCase(directory, 
										  t_start=args.times[0], 
										  t_end=args.times[1])

	# read coefficients and compute mean values of cuIBM simulation
	if args.cuibm_path:
		cases['cuibm'] = CuIBMCase(cases['cuibm'], 
								   t_start=args.times[0], t_end=args.times[1])

	# read drag coefficient from Koumoutsakos and Leonard (1995)
	if args.kl1995:
		# path of data file
		kl1995_path = ('/home/mesnardo/flying_snake_openfoam/resources/'
						'cylinder_drag_coefficient_Re550_'
						'koumoutsakos_leonard_1995.dat')
		cases['kl1995'] = Case(kl1995_path)
		with open(cases['kl1995'].path, 'r') as infile:
			cases['kl1995'].t, cases['kl1995'].cd = numpy.loadtxt(infile, 
																dtype=float,
																delimiter='\t',
																unpack=True)
			cases['kl1995'].t *= 0.5	# to use the same time-scale

	# write mean force coefficients and Strouhal numbers in log file
	if args.save:
		logging.info('[case] %s' % cases['main'].path)
		cd_mean, cl_mean = cases['main'].cd_mean, cases['main'].cl_mean
		cd_min = cases['main'].cd_mean - cases['main'].cd_min
		cd_max = cases['main'].cd_max - cases['main'].cd_mean
		cl_min = cases['main'].cl_mean - cases['main'].cl_min
		cl_max = cases['main'].cl_max - cases['main'].cl_mean
		logging.info('\tcd = %f (-%f, +%f)' % (cd_mean, cd_min, cd_max))
		logging.info('\tcl = %f (-%f, +%f)' % (cl_mean, cl_min, cl_max))
		logging.info('\tSt = %f' % cases['main'].strouhal)
		for case in cases['others']:
			logging.info('[case] %s' % case.path)
			cd_mean, cl_mean = case.cd_mean, case.cl_mean
			cd_min = case.cd_mean - case.cd_min
			cd_max = case.cd_max - case.cd_mean
			cl_min = case.cl_mean - case.cl_min
			cl_max = case.cl_max - case.cl_mean
			logging.info('\tcd = %f (-%f, +%f)' % (cd_mean, cd_min, cd_max))
			logging.info('\tcl = %f (-%f, +%f)' % (cl_mean, cl_min, cl_max))
			logging.info('\tSt = %f' % case.strouhal)

	# plot force coefficients
	plot_coefficients(cases, args)


if __name__ == '__main__':
	main()
