#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_force_coefficients.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot aerodynamic coefficients of an OpenFoam case


import os
import sys
import argparse
import collections

import numpy
import pandas
from scipy import signal
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the lift and drag '
									'coefficients for a given simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--case', dest='directory', type=str, 
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
	parser.add_argument('--order', dest='order', type=int, default=5,
						help='number of neighbors to compare with to get limit '
							 'indices')
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
	parser.add_argument('--compare', dest='other_directories', type=str, 
						nargs='+',default=[],
						help='directories of other cases for comparison')
	parser.add_argument('--cuibm', dest='cuibm', type=str, default=None,
						help='path of cuIBM force coefficients for comparison')
	parser.add_argument('--reference', dest='reference', type=str, default=None,
						help='case reference to compute relative errors')
	parser.add_argument('--kl1995', dest='kl1995', action='store_true',
						help='plots instantaneous drag coefficient from '
							 'Koumoutsakos and Leonard (1995)')
	parser.set_defaults(save=True, lift=True, drag=True)
	# parse the command-line
	return parser.parse_args()


def relative_change(current, reference):
	"""Computes the relative change as a percentage.
	
	Arguments
	---------
	current -- current value.
	reference -- value used as a reference.
	"""
	return 100.0 * (current-reference) / abs(reference)


class ForceCoefficient(object):
	"""Contains the values of a force coefficient, 
	as well as its means and fluctuations and errors relative to a reference.
	"""
	def __init__(self, name):
		"""Initializes the force coefficient and its name.
		
		Arguments
		---------
		name -- name of the force coefficient.
		"""
		self.name = name
		self.values = numpy.empty(0)

	def get_mean(self, i_start, i_end):
		"""Computes the mean force coefficient.
		
		Arguments
		---------
		i_start, i_end -- indices of the interval to work with.
		"""
		self.mean = numpy.mean(self.values[i_start:i_end])

	def get_fluctuations(self, i_start, i_end):
		"""Computes the fluctuations about the mean.
		
		Arguments
		---------
		i_start, i_end --indices of the interval to work with.
		"""
		minimum = self.values[i_start:i_end].min()
		maximum = self.values[i_start:i_end].max()
		self.fluctuations = [minimum, maximum] - self.mean

	def get_relative_errors(self, reference):
		"""Computes the relative errors.
		
		Arguments
		---------
		reference -- force coefficient to use as a reference.
		"""
		self.errors = [relative_change(self.mean, reference.mean),
					   relative_change(self.fluctuations[0], 
					   				   reference.fluctuations[0]),
					   relative_change(self.fluctuations[1], 
					   				   reference.fluctuations[1])]


class Case(object):
	"""Contains information about the simulation."""
	def __init__(self, path, args):
		"""Stores the path and the legend of the simulation,
		as well as information about the interval to work with.
		
		Arguments
		---------
		path -- path of the simulation.
		args -- command-line arguments.
		"""
		self.path = path
		self.legend = os.path.basename(self.path)
		self.t_start = (None if not args.times[0] else args.times[0])
		self.t_end = (None if not args.times[1] else args.times[1])
		self.order = (None if not args.order else args.order)

	def get_limit_indices(self):
		"""Computes the indices in the array of the time-limits."""
		if self.t_start and self.t_end:
			self.i_start = numpy.where(self.t >= self.t_start)[0][0]
			self.i_end = numpy.where(self.t >= self.t_end)[0][0]-1
		else:
			minima = signal.argrelextrema(self.cl.values, numpy.less_equal, 
										  order=self.order)[0][:-1]
			minima = minima[numpy.append(True, 
										 (minima[1:]-minima[:-1])>self.order)]
			self.i_start, self.i_end = minima[-2], minima[-1]

	def get_mean_coefficients(self):
		"""Computes the mean force coefficients."""
		self.cd.get_mean(self.i_start, self.i_end)
		self.cl.get_mean(self.i_start, self.i_end)

	def get_fluctuations(self):
		"""Computes the fluctuations of the force coefficients."""
		self.cd.get_fluctuations(self.i_start, self.i_end)
		self.cl.get_fluctuations(self.i_start, self.i_end)

	def get_relative_errors(self, reference):
		"""Computes the relative error of the mean and fluctuations 
		of each force coefficient.
		
		Arguments
		---------
		reference -- case reference to compute relative errors.
		"""
		self.cd.get_relative_errors(reference.cd)
		self.cl.get_relative_errors(reference.cl)

	def get_strouhal_number(self, is_strouhal):
		"""Computes the Strouhal number.
		
		Arguments
		---------
		is_strouhal -- boolean to compute or not the Strouhal number.
		"""
		self.strouhal = 0.0
		if is_strouhal:
			self.strouhal = 1.0/(self.t[self.i_end] - self.t[self.i_start])


class OpenFoamCase(Case):
	"""Contains force coefficients of an OpenFoam simulation."""
	def __init__(self, directory, args):
		"""Reads the force coefficients and computes the mean coefficients, 
		the fluctuations and the Strouhal number.
		
		Arguments
		---------
		directory -- directory of the simulation.
		args -- command-line arguments.
		"""
		Case.__init__(self, directory, args)
		# read force coefficients
		self.read_coefficients()
		# compute time-interval indices
		self.get_limit_indices()
		# compute mean coefficients
		self.get_mean_coefficients()
		# compute fluctuations
		self.get_fluctuations()
		# compute Strouhal number
		self.get_strouhal_number((True if self.t_start == None else False))

	def read_coefficients(self):
		"""Reads force coefficients from files."""
		# initialization
		self.t = numpy.empty(0)
		self.cd, self.cl = ForceCoefficient('cd'), ForceCoefficient('cl')
		# read force coefficients from files and append to arrays
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
			self.cd.values = numpy.append(self.cd.values, cd_tmp)
			self.cl.values = numpy.append(self.cl.values, cl_tmp)


class CuIBMCase(Case):
	"""Contains force coefficients of a cuIBM  simulation."""
	def __init__(self, path, args):
		"""Reads the force coefficients and computes the mean coefficients, 
		the fluctuations and the Strouhal number.
		
		Arguments
		---------
		path -- path of the cuIBM force coefficients file.
		args -- command-line arguments.
		"""
		Case.__init__(self, path, args)
		self.legend = 'cuIBM'
		# read force coefficients
		self.read_coefficients()
		# compute time-interval indices
		self.get_limit_indices()
		# compute mean coefficients
		self.get_mean_coefficients()
		# compute fluctuations
		self.get_fluctuations()
		# compute Strouhal number
		self.get_strouhal_number((True if self.t_start == None else False))

	def read_coefficients(self):
		"""Reads force coefficients from cuIBM results."""
		# read the file
		self.cd, self.cl = ForceCoefficient('cd'), ForceCoefficient('cl')
		with open(self.path, 'r') as infile:
			self.t, self.cd.values, self.cl.values = numpy.loadtxt(infile, 
																dtype=float, 
													 			delimiter='\t', 
													 			unpack=True)
		# cuIBM script does not include the coefficient 2.0
		self.cd.values *= 2.0
		self.cl.values *= 2.0


class KL1995Case(Case):
	"""Contain drag coefficient results from Koumoutsakos and Leonard (1995)."""
	def __init__(self, path):
		"""Get instantaneous drag coefficient 
		reported by Koumoutsakos and Leonard (1995).
		"""
		Case.__init__(self, path, args)
		self.read_drag_coefficients()

	def read_drag_coefficients(self):
		"""Reads and stores the intantaneous drag coefficient."""
		with open(self.path, 'r') as infile:
			self.cd = ForceCoefficient('cd')
			self.t, self.cd.values = numpy.loadtxt(infile, dtype=float, 
												   delimiter='\t', unpack=True)
			self.t *= 0.5	# to use the same time-scale


def plot_coefficients(cases, args):
	"""Plots force coefficients from different simulations.
	
	Arguments
	--------
	cases -- dictionary that contains info about all simulations to plot.
	args -- command-line arguments
	"""
	# figure parameters
	pyplot.figure(figsize=(10, 6))
	pyplot.grid(True)
	pyplot.xlabel('time', fontsize=16)
	pyplot.ylabel('force coefficients', fontsize=16)
	# plot the main OpenFoam force coefficients
	if args.drag:
		pyplot.plot(cases[args.directory].t, cases[args.directory].cd.values, 
					label=r'$C_d$ - %s' % cases[args.directory].legend, 
					color='r', ls='-', lw=2)
	if args.lift:
		pyplot.plot(cases[args.directory].t, cases[args.directory].cl.values, 
					label=r'$C_l$ - %s' % cases[args.directory].legend, 
					color='b', ls='-', lw=2)
	# plot other OpenFoam force coefficients
	colors = ['g', 'c', 'm', 'y']
	for i, directory in enumerate(args.other_directories):
		if args.drag:
			pyplot.plot(cases[directory].t, cases[directory].cd.values, 
						label=r'$C_d$ - %s' % cases[directory].legend,
						color=colors[i], ls='-', lw=1)
		if args.lift:
			pyplot.plot(cases[directory].t, cases[directory].cl.values,
						label=r'$C_l$ - %s' % cases[directory].legend,
						color=colors[i], ls='--', lw=1)
	# plot cuIBM force coefficients
	if cases[args.cuibm]:
		if args.drag:
			pyplot.plot(cases[args.cuibm].t, cases[args.cuibm].cd.values,
						label=r'$C_d$ - cuIBM', color='k', ls='-', lw=1)
		if args.drag:
			pyplot.plot(cases[args.cuibm].t, cases[args.cuibm].cl.values,
						label=r'$C_l$ - cuIBM', color='k', ls='--', lw=1)
	# plot Koumoutsakos and Leonard (1995) drag coefficients
	if args.kl1995:
		pyplot.plot(cases['kl1995'].t, cases['kl1995'].cd,
					label=r'$C_d$ - Koumoutsakos and Leonard (1995)', 
					color='k', lw=0, marker='o', markersize=6)
	# define limits of the figure
	x_min = (cases[args.directory].t[0] if not args.limits else args.limits[0])
	x_max = (cases[args.directory].t[-1] if not args.limits else args.limits[1])
	y_min = (-2.0 if not args.limits else args.limits[2])
	y_max = (+2.0 if not args.limits else args.limits[3])
	pyplot.xlim(x_min, x_max)
	pyplot.ylim(y_min, y_max)
	# add legend to the figure
	pyplot.legend(loc='upper left', prop={'size': 12},
				  bbox_to_anchor=(0.0, -0.05))
	# save the figure as a .PNG file
	if args.save:
		# create images folder if not existing
		images_directory = '%s/images' % args.directory
		if not os.path.isdir(images_directory):
			os.makedirs(images_directory)
		pyplot.savefig('%s/%s.png' % (images_directory, args.image_name),
					   bbox_inches='tight')
	# display figure
	if args.show:
		pyplot.show()


def print_coefficients(cases, args):
	"""Prints the mean, fluctuations and relative errors 
	of the different simulations.

	Arguments
	---------
	cases -- dictionary containing the different cases.
	args -- command-line arguments.
	"""
	# get the legend for each simulation
	if args.legend:
		cases[args.directory].legend = args.legend[0]
		for i, directory in enumerate(args.other_directories):
			cases[directory].legend = args.legend[i+1]
	# compute relative errors
	for key, case in cases.iteritems():
		if args.reference:
			case.get_relative_errors(cases[args.reference])

	# write data in csv file
	csv_path = '%s/postProcessing/%s.csv' % (args.directory, args.image_name)
	with open(csv_path, 'w') as outfile:
		# write drag coefficients
		outfile.write('simulation\tcd\tmin\tmax\n')
		for key, case in cases.iteritems():
			if case.path == args.reference:
				outfile.write('%s\t%.4f\t%.4f\t%.4f\n' % (case.legend,
												   	case.cd.mean,
												   	case.cd.fluctuations[0],
												   	case.cd.fluctuations[1]))
			else:
				case.get_relative_errors(cases[args.reference])
				outfile.write('%s\t%.4f (%+.1f%%)\t'
							  '%.4f (%+.1f%%)\t%.4f (%+.1f%%)\n'
							  % (case.legend, 
							 	 case.cd.mean, case.cd.errors[0],
								 case.cd.fluctuations[0], case.cd.errors[1],
								 case.cd.fluctuations[1], case.cd.errors[2]))
		# write lift coefficients
		outfile.write('simulation\tcl\tmin\tmax\n')
		for key, case in cases.iteritems():
			if case.path == args.reference:
				outfile.write('%s\t%.4f\t%.4f\t%.4f\n' % (case.legend,
												   	case.cl.mean,
												   	case.cl.fluctuations[0],
												   	case.cl.fluctuations[1]))
			else:
				case.get_relative_errors(cases[args.reference])
				outfile.write('%s\t%.4f (%+.1f%%)\t'
							  '%.4f (%+.1f%%)\t%.4f (%+.1f%%)\n'
							  % (case.legend, 
							 	 case.cl.mean, case.cl.errors[0],
								 case.cl.fluctuations[0], case.cl.errors[1],
								 case.cl.fluctuations[1], case.cl.errors[2]))
	
	# create dataframes for force coefficients using pandas
	df_cd = pandas.read_csv(csv_path, 
							delimiter='\t', comment='#', nrows=len(cases),
							header=0, index_col=0)
	df_cl = pandas.read_csv(csv_path, 
							delimiter='\t', comment='#', nrows=len(cases),
							header=len(cases)+1, index_col=0)
	if args.times[0]:
		print ('\nAveraging the instantaneous force coefficients\n'
			   'between %g and %g seconds of flow simulation:' 
			   % (args.times[0], args.times[1]))
		print '\n', df_cd, '\n\n', df_cl, '\n'
	else:
		# write results in csv file
		with open(csv_path, 'a') as outfile:
			outfile.write('simulation\tSt\n')
			for key, case in cases.iteritems():
				if case.path == args.reference:
					outfile.write('%s\t%.4f\n' % (case.legend, case.strouhal))
				else:
					# compute relative error in Strouhal number
					error_strouhal = relative_change(case.strouhal, 
												cases[args.reference].strouhal)
					outfile.write('%s\t%.4f (%+.1f%%)\n' % (case.legend, 
														case.strouhal, 
														error_strouhal))
		# create dataframe for Strouhal number using pandas
		df_st = pandas.read_csv(csv_path,
								delimiter='\t', comment='#', nrows=len(cases),
								header=2*(len(cases)+1), index_col=0)
		# display dataframes
		print ('\nAveraging the instantaneous force coefficients\n'
			   'over the last period:')
		print '\n', df_cd, '\n\n', df_cl, '\n\n', df_st, '\n'

	
def main():
	"""Plots aerodynamic coefficients."""
	# parse the command-line
	args = read_inputs()

	# dictionary containing the force coefficients of each case
	cases = collections.OrderedDict()

	# read and analyze force coefficients from cuIBM results
	if args.cuibm:
		cases[args.cuibm] = CuIBMCase(args.cuibm, args)

	# read and analyze force coefficients from main OpenFoam case
	cases[args.directory] = OpenFoamCase(args.directory, args)

	# read and analyze force coefficients from other OpenFoam cases
	for directory in args.other_directories:
		cases[directory] = OpenFoamCase(directory, args)

	# read drag coefficient from Koumoutsakos and Leonard (1995)
	if args.kl1995:
		# path of data file
		kl1995_path = ('/home/mesnardo/flying_snake_openfoam/resources/'
						'cylinder_drag_coefficient_Re550_'
						'koumoutsakos_leonard_1995.dat')
		cases['kl1995'] = KL1995Case(kl1995_path)

	# print force coefficients
	print_coefficients(cases, args)

	# plot force coefficients
	plot_coefficients(cases, args)

if __name__ == '__main__':
	main()
