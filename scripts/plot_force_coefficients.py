#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_force_coefficients.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot aerodynamic coefficients of an OpenFoam case


import os
import argparse
import datetime

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
							 ' and the Strouhal number')
	parser.add_argument('--end', dest='end', type=float, default=None,
						help='ending-time to compute mean coefficients '
							 ' and the Strouhal number')
	parser.add_argument('--limits', dest='limits', type=float, nargs='+',
						default=[None, None, None, None],
						help='time-limits to plot force coefficients')
	parser.add_argument('--cuibm', dest='cuibm_path', type=str, default=None,
						help='path of cuIBM force coefficients for comparison')
	parser.add_argument('--kl1995', dest='kl1995', action='store_true',
						help='plots instantaneous drag coefficient from '
							 'Koumoutsakos and Leonard (1995)')
	parser.add_argument('--name', dest='image_name', type=str, 
						default='force_coefficients',
						help='name of the file generated (no extension)')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure as a .png file')
	parser.add_argument('--compare', dest='other_cases', type=str, nargs='+',
						help='directories of other cases for comparison')
	parser.add_argument('--legend', dest='legend', type=str, nargs='+',
						help='legend for each simulation to plot')
	parser.add_argument('--no-lift', dest='lift', action='store_false',
						help='does not plot the lift coefficients')
	parser.add_argument('--no-drag', dest='drag', action='store_false',
						help='does not plot the drag coefficients')
	parser.set_defaults(save=True, lift=True, drag=True)
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

	def get_strouhal_number(self):
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
		self.get_strouhal_number()

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
	print args.limits
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

	# write list of command-line arguments in a log file
	log_path = ('%s/%s.log' % (args.case_directory, 
		  					   os.path.splitext(os.path.basename(__file__))[0]))
	if args.save:
		with open(log_path, 'w') as outfile:
			outfile.write('%s\n%s\n' % (str(datetime.datetime.now()), 
										str(args)))

	# store case directories
	cases = {'main': args.case_directory,
			 'others': ([] if not args.other_cases else args.other_cases),
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

	# read drag coefficient from Koumoutsakos and Leonard (1995)
	if args.kl1995:
		kl1995_path = ('$FLYING_SNAKE_OPENFOAM/resources/'
						'cylinder_drag_coefficient_Re550_'
						'koumoutsakos_leonard_1995.dat')
		kl1995_path = ('/home/mesnardo/flying_snake_openfoam/resources/'
						'cylinder_drag_coefficient_Re550_'
						'koumoutsakos_leonard_1995.dat')
		cases['kl1995'] = Case(kl1995_path)
		with open(cases['kl1995'].path, 'r') as infile:
			cases['kl1995'].t, cases['kl1995'].cd = numpy.loadtxt(infile, 
																dtype=float,
																delimiter='\t',
																unpack=True)
			cases['kl1995'].t *= 0.5

	# write mean force coefficients and Strouhal numbers in log file
	if args.save:
		with open(log_path, 'a') as outfile:
			outfile.write('\ncase: %s\n' % cases['main'].path)
			outfile.write('cd = %f\n' % cases['main'].cd_mean)
			outfile.write('cl = %f\n' % cases['main'].cl_mean)
			outfile.write('St = %f\n' % cases['main'].strouhal)
			for case in cases['others']:
				outfile.write('\ncase: %s\n' % case.path)
				outfile.write('cd = %f\n' % case.cd_mean)
				outfile.write('cl = %f\n' % case.cl_mean)
				outfile.write('St = %f\n' % case.strouhal)


	# plot force coefficients
	plot_coefficients(cases, args)


if __name__ == '__main__':
	main()
