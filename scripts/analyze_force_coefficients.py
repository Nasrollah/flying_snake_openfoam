#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/wip_force_coefficients_convergence.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Evaluates the convergence of steady or pseudo-steady simulations.


import os
import argparse


import numpy
from scipy import signal
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create parser
	parser = argparse.ArgumentParser(description='Evaluates convergence of '
						'force coefficients of steady and pseudo-steady '
						'OPenFoam/icoFoam simulations',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill parser wih arguments
	parser.add_argument('--case', dest='case_directory', type=str,
						default=os.getcwd(),
						help='directory of the OpenFoam case')
	parser.add_argument('--limits-cd', dest='limits_cd', type=float, nargs='+',
						default=[None, None, 1.36, 1.39],
						help='limits to plot the drag coefficient')
	parser.add_argument('--limits-cl', dest='limits_cl', type=float, nargs='+',
						default=[None, None, -0.40, 0.40],
						help='limits to plot the lift coefficient')
	parser.add_argument('--no-frequency', dest='plot_frequency', 
						action='store_false',
						help='does not plot the frequency of the force')
	parser.add_argument('--no-mean', dest='plot_mean', action='store_false',
						help='does not plot the means')
	parser.add_argument('--no-extrema', dest='plot_extrema', 
						action='store_false',
						help='does not plot the extrema')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figures')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the figures')
	parser.set_defaults(plot_frequency=True, plot_mean=True, plot_extrema=True, 
						save=True)
	return parser.parse_args()


def read_force_coefficients(case_directory):
	"""Reads the instantaneous force coefficients from a given OpenFoam 
	simulation directory.
	
	Arguments
	---------
	case_directory -- directory of the simulation folder.

	Returns
	-------
	t -- the discrete time values.
	cd, cl -- discrete drag and lift coefficient values.
	"""
	t, cd, cl = numpy.empty(0), numpy.empty(0), numpy.empty(0)
	forces_directory = '%s/postProcessing/forceCoeffs' % case_directory
	folders = sorted(os.listdir(forces_directory))
	for folder in folders:
		force_path = '%s/%s/forceCoeffs.dat' % (forces_directory, folder)
		with open(force_path, 'r') as infile:
			t_tmp, cd_tmp, cl_tmp = numpy.loadtxt(infile, dtype=float, 
												  usecols=(0, 2, 3), 
												  delimiter='\t', 
												  unpack=True)
		t = numpy.append(t, t_tmp)
		cd = numpy.append(cd, cd_tmp)
		cl = numpy.append(cl, cl_tmp)
	return t, cd, cl


class ForceCoefficient:
	"""Contains info about a force coefficient."""
	def __init__(self, name, t, f):
		"""Stores name and instantaneous force coefficient.
		
		Arguments
		---------
		name -- name of the force coefficient.
		t -- the time.
		f -- the instantaneous force coefficient.
		"""
		self.name = name
		self.t = t
		self.f = f

	def get_extrema_indices(self):
		"""Computes extrema indices of the force coefficient."""
		self.minima = signal.argrelextrema(self.f, 
										   numpy.less_equal, order=5)[0][:-1]
		self.maxima = signal.argrelextrema(self.f, 
										   numpy.greater_equal, order=5)[0][:-1]
		# remove close values
		self.minima = numpy.delete(self.minima, 
								   [i+1 for i in xrange(self.minima.size-1) 
								   	if (self.minima[i+1] == self.minima[i]+1)])
		self.maxima = numpy.delete(self.maxima, 
								   [i+1 for i in xrange(self.maxima.size-1) 
								   	if (self.maxima[i+1] == self.maxima[i]+1)])

	def get_periods(self):
		"""Calculates the periods over the time of the simulation 
		using the minima, and computes the Strouhal number reached.
		"""
		self.periods = numpy.empty(self.minima.size-2, dtype=float)
		for i in xrange(self.periods.size):
			self.periods[i] = self.t[self.minima[i+1]] - self.t[self.minima[i]]
		self.last_period = self.periods[-1]
		self.strouhal = 1.0 / self.last_period

	def get_means(self):
		"""Computes the mean coefficient over each period using the minima,
		and stores the last mean coefficient reached.
		"""
		self.means = numpy.empty(self.minima.size-2, dtype=float)
		for i in xrange(self.means.size):
			slice_f = self.f[self.minima[i]:self.minima[i+1]]
			self.means[i] = (self.f[self.minima[i]:self.minima[i+1]].sum() 
							 / self.f[self.minima[i]:self.minima[i+1]].size)
		self.last_mean = self.means[-1]

	def print_info(self):
		"""Displays info about the force into the terminal."""
		print '\nname: %s' % self.name
		print 'mean: %.6f' % self.last_mean
		print 'Strouhal: %.6f\n' % self.strouhal

	def plot(self, args, limits=[None, None, None, None]):
		"""Plots the instantaneous force coefficients and its peaks.
		
		Arguments
		---------
		args -- command-line arguments.
		limits -- limits of the plot (default: [None, None, None, None]).
		"""
		fig, ax1 = pyplot.subplots()
		ax1.set_xlabel(r'$t$', fontsize=18)
		ax1.set_ylabel(self.name, fontsize=16)
		ax1.plot(self.t, self.f, label=self.name, color='k', lw=1, ls='-')
		# means
		if args.plot_mean:
			t_half = 0.5 * (self.t[self.minima[:-2]]+self.t[self.minima[1:-1]])
			ax1.plot(t_half, self.means, label='mean', 
				 	 color='r', lw=0, marker='o', markersize=6)
			ax1.axhline(self.last_mean, label=None, 
						color='r', lw=1, ls='--')
		# extrema
		if args.plot_extrema:
			# minima
			ax1.plot(self.t[self.minima], self.f[self.minima], label='minima',
					 color='g', lw=0, marker='o', markersize=6)
			ax1.axhline(self.f[self.minima[-1]], label=None,
						color='g', lw=1, ls='--')
			# maxima
			ax1.plot(self.t[self.maxima], self.f[self.maxima], label='maxima',
					 color='b', lw=0, marker='o', markersize=6)
			ax1.axhline(self.f[self.maxima[-1]], label=None,
						color='b', lw=1, ls='--')
		ax1.set_xlim(limits[0],limits[1])
		ax1.set_ylim(limits[2], limits[3])
		ax1.legend(loc='upper left', prop={'size': 12})
		# periods
		if args.plot_frequency:
			ax2 = ax1.twinx()
			ax2.set_ylabel('Strouhal', color='m', fontsize=16)
			ax2.plot(t_half, 1.0/self.periods, 
					 color='m', lw=0, marker='o', markersize=6)
			ax2.axhline(self.last_period, color='m', lw=1, ls='--')
			ax2.set_ylim(0.9*self.strouhal, 1.1*self.strouhal)
		# save the plot
		if args.save:
			# create images folder if not existing
			images_directory = '%s/images' % args.case_directory
			if not os.path.isdir(images_directory):
				os.makedirs(images_directory)
			pyplot.savefig('%s/analyze%s.png' % (images_directory, self.name))

	def analyze(self, args, limits=[None, None, None, None]):
		"""Computes means, periods, Strouhal number and plots the resuls in a 
		figure.
		
		Arguments
		---------
		args -- command-line arguments
		limits -- limits of the plot (default: [None, None, None, None]).
		"""
		self.get_extrema_indices()
		self.get_periods()
		self.get_means()
		self.print_info()
		self.plot(args, limits=limits)


def main():
	"""Evaluates the convergence of steady or pseudo-steady simulations."""
	# parse command-line
	args = read_inputs()

	print '\n[case] %s' % args.case_directory

	# read force coefficients
	t, cd, cl = read_force_coefficients(args.case_directory)

	# create force coefficient objects
	cd = ForceCoefficient('Cd', t, cd)
	cl = ForceCoefficient('Cl', t, cl)

	# analyze the force coefficients
	cd.analyze(args, limits=args.limits_cd)
	cl.analyze(args, limits=args.limits_cl)

	# display figures
	if args.show:
		pyplot.show()


if __name__ == '__main__':
	main()
