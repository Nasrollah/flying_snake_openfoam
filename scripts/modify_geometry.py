#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/modify_geometry.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Discretization, rotation, translation and scaling 
#			   of a given geometry


import os
import argparse
import math

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Modify a given geometry',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--coordinates', dest='coordinates_path', type=str,
						help='path of the coordinates file')
	parser.add_argument('-n', dest='n', type=int, default=None,
						help='number of points for discretization')
	parser.add_argument('--cl', dest='cl', type=float, default=None,
						help='characteristic-length on the geometry')
	parser.add_argument('--scale', '-s', dest='scale', type=float, default=1.0,
						help='scaling factor of the new geometry')
	parser.add_argument('--rotation', '-r', dest='rotation', type=float, 
						default=0.0,
						help='angle of rotation in degrees')
	parser.add_argument('--center', '-c', dest='center', type=float,
						default=[None, None],
						help='center of rotation')
	parser.add_argument('--translation', '-t', dest='translation', type=float, 
						nargs='+', default=[0.0, 0.0, 0.0],
						help='displacement in the x- and y- directions')
	parser.add_argument('--name', dest='name', type=str, default='new_body',
						help='name of the output file (without extension)')
	parser.add_argument('--extension', dest='extension', type=str, 
						default='bdy',
						help='extension of the file generated')
	parser.add_argument('--save-dir', dest='save_dir', type=str, 
						default=os.getcwd(),
						help='directory where to save the new coordinates file')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the initial and current geometries')
	return parser.parse_args()


class Geometry:
	"""Definition of the geometry."""
	def __init__(self, file_path):
		"""Initializes the geometry by reading the coordinates file.
		
		Arguments
		---------
		file_name -- path of the coordinates file.
		"""
		self.read(file_path)
		self.get_center_mass()

	def read(self, file_path):
		"""Stores the coordinates of the geometry into arrays.
		
		Arguments
		--------
		file_path -- path of the coordinates file.
		"""
		with open(file_path, 'r') as infile:
			self.x, self.y = numpy.loadtxt(infile, dtype=float, 
						  	  			   delimiter='\t', unpack=True, 
										   skiprows=1)
		# store initial coordinates
		self.x_old, self.y_old = self.x.copy(), self.y.copy()

	def write(self, file_path):
		"""Writes the coordinates of the geometry into a file.
		
		Arguments
		---------
		file_path -- path of the output file.
		"""
		with open(file_path, 'w') as outfile:
			outfile.write('%d\n' % self.x.size)
			numpy.savetxt(outfile, numpy.c_[self.x, self.y],
					   	  fmt='%.6f', delimiter='\t')

	def get_distance(self, x_start, y_start, x_end, y_end):
		"""Returns the distance between two points.
		
		Arguments
		---------
		x_start, y_end -- coordinates of the first point.
		x_end, y_end -- coordinates of the second point.
		"""
		return math.sqrt((x_end-x_start)**2 + (y_end-y_start)**2)

	def get_perimeter(self):
		"""Returns the perimeter of the geometry."""
		x, y = numpy.append(self.x, self.x[0]), numpy.append(self.y, self.y[0])
		return numpy.sum(numpy.sqrt((x[1:]-x[0:-1])**2+(y[1:]-y[0:-1])**2))

	def get_center_mass(self):
		"""Computes the center of mass."""
		self.x_cm = self.x.sum()/self.x.size
		self.y_cm = self.y.sum()/self.y.size

	def rotation(self, aoa=0.0, x_rot=None, y_rot=None):
		"""Rotates the geometry.
		
		Arguments
		---------
		aoa -- angle of rotation in degrees (default 0.0).
		x_rot, y_rot -- location of the center of rotation (default None, None).
		"""
		if not (x_rot and y_rot):
			x_rot, y_rot = self.x_cm, self.y_cm
		aoa *= math.pi/180.
		x_tmp = x_rot + (self.x-x_rot)*math.cos(aoa) \
					  - (self.y-y_rot)*math.sin(aoa)
		y_tmp = y_rot + (self.x-x_rot)*math.sin(aoa) \
					  + (self.y-y_rot)*math.cos(aoa)
		self.x, self.y = x_tmp, y_tmp
		self.get_center_mass()

	def translation(self, x_trans=0.0, y_trans=0.0):
		"""Translates the geometry.
		
		Arguments
		---------
		x_trans, y_trans -- x- and y- displacements (default 0.0, 0.0).
		"""
		self.x += x_trans
		self.y += y_trans
		self.get_center_mass()

	def scale(self, ratio=1.0):
		"""Scales the geometry.
		
		Arguments
		---------
		ratio -- scaling ratio (default 1.0).
		"""
		self.x = self.x_cm + ratio*(self.x - self.x_cm)
		self.y = self.y_cm + ratio*(self.y - self.y_cm)
	
	def interpolation(self, x_start, y_start, x_end, y_end, cl):
		"""Computes the coordinates of a point 
		by interpolation between two given points given a distance.
		
		Arguments
		---------
		x_start, y_start -- coordinates of the starting point.
		x_end, y_end -- coordinates of the ending point.
		cl -- length between the starting point and the interpolated one.

		Returns
		-------
		x_target, y_target -- coordinates of the interpolated point.
		"""
		length = self.get_distance(x_start, y_start, x_end, y_end)
		x_target = x_start + cl/length*(x_end-x_start)
		y_target = y_start + cl/length*(y_end-y_start)
		return x_target, y_target

	def projection(self, x_start, y_start, x_tmp, y_tmp, x_end, y_end, cl):
		"""Computes the coordinates of a point
		by projection onto the segment [(x_tmp, y_tmp), (x_end, y_end)]
		such that the distance between (x_start, y_start) and the new point
		is the given cahracteristic length.

		Arguments
		---------
		x_start, y_start -- coordinates of the starting point.
		x_tmp, y_tmp -- coordinates of the intermediate point.
		x_end, y_end -- coordinates of the ending point.
		cl -- characteristic length.

		Returns
		-------
		x_target, y_target -- coordinates of the projected point.
		"""
		tol = 1.0E-06
		if abs(y_end-y_tmp) >= tol:
			# solve for y
			# coefficients of the second-order polynomial
			a = (x_end-x_tmp)**2 + (y_end-y_tmp)**2
			b = 2.0*( (x_end-x_tmp)*( y_tmp*(x_start-x_end) 
									+ y_end*(x_tmp-x_start) ) 
					- y_start*(y_end-y_tmp)**2 )
			c = (y_start**2-cl**2)*(y_end-y_tmp)**2 \
				+ (y_tmp*(x_start-x_end) + y_end*(x_tmp-x_start))**2
			# solve the second-order polynomial: ay^2 + by + c = 0
			y = numpy.roots([a, b, c])
			# test if the point belongs to the segment
			test = (y_tmp <= y[0] <= y_end or y_end <= y[0] <= y_tmp)
			y_target = (y[0] if test else y[1])
			x_target = x_tmp + (x_end-x_tmp)/(y_end-y_tmp)*(y_target-y_tmp)
		else:
			# solve for x
			# coefficients of the second-order polynomial
			a = (x_end-x_tmp)**2 + (y_end-y_tmp)**2
			b = 2.0*( (x_end-x_tmp)*(y_tmp-y_start)*(y_end-y_tmp) 
					 - x_start*(x_end-x_tmp)**2 
					 - x_tmp*(x_end-x_tmp)**2 )
			c = (x_end-x_tmp)**2*((y_tmp-y_start)**2+x_start**2-cl**2) \
				+ x_tmp**2*(y_end-y_tmp)**2 \
				- 2*x_tmp*(x_end-x_tmp)*(y_tmp-y_start)*(y_end-y_tmp)
			# solve the second-order polynomial: ax^2 + bx + c = 0
			x = numpy.roots([a, b, c])
			# test if the point belongs to the segment
			test = (x_tmp <= x[0] <= x_end or x_end <= x[0] <= x_tmp)
			x_target = (x[0] if test else x[1])
			y_target = y_tmp + (y_end-y_tmp)/(x_end-x_tmp)*(x_target-x_tmp)
		return x_target, y_target

	def discretize(self, n=None, cl=None):
		"""Discretizes the geometry 
		given a characteristic length or a number of points.
		
		Arguments
		---------
		n -- number of points (default None).
		cl -- characteristic length (default None).
		"""
		# calculate either the numebr of points or characteristic length
		if n and not cl:
			cl = self.get_perimeter()/n
		elif cl and not n:
			n = int(self.get_perimeter()/cl)
		elif not (n and cl):
			return

		# exit function if same discretization
		if n == self.x.size:
			return

		# copy coordinates and initialize new ones
		x_old = numpy.append(self.x, self.x[0]) 
		y_old = numpy.append(self.y, self.y[0])
		x_new, y_new = numpy.empty(n, dtype=float), numpy.empty(n, dtype=float)
		# first element
		x_new[0], y_new[0] = x_old[0], y_old[0]

		I = 0
		tol = 1.0E-06    # tolerance for interpolation
		for i in xrange(n-1):
			x_start, y_start = x_new[i], y_new[i]
			x_end, y_end = x_old[I+1], y_old[I+1]
			distance = self.get_distance(x_start, y_start, x_end, y_end)
			if cl-distance <= tol:
				# interpolation method
				x_new[i+1], y_new[i+1] = self.interpolation(x_start, y_start, 
															x_end, y_end, 
															cl)
			else:
				# projection method
				while I < x_old.size-2 and cl-distance > tol:
					I += 1
					x_tmp, y_tmp = x_end, y_end
					x_end, y_end = x_old[I+1], y_old[I+1]
					distance = self.get_distance(x_start, y_start, x_end, y_end)
				x_new[i+1], y_new[i+1] = self.projection(x_start, y_start, 
														 x_tmp, y_tmp, 
														 x_end, y_end, 
														 cl)
		# store the new discretization
		self.x, self.y = x_new.copy(), y_new.copy()

	def plot(self):
		"""Plots the geometry."""
		pyplot.figure()
		pyplot.grid(True)
		pyplot.xlabel(r'$x$', fontsize=18)
		pyplot.ylabel(r'$y$', fontsize=18)
		pyplot.plot(numpy.append(self.x_old, self.x_old[0]), 
				 	numpy.append(self.y_old, self.y_old[0]),
				 	label='initial',
				 	color='k', ls='-', lw=2, marker='o', markersize=6)
		pyplot.plot(numpy.append(self.x, self.x[0]), 
				 	numpy.append(self.y, self.y[0]),
				 	label='current',
				 	color='r', ls='-', lw=2, marker='o', markersize=6)
		pyplot.axis('equal')
		pyplot.legend(loc='best', prop={'size': 16})
		pyplot.show()


def main():
	"""Applies discretization, scaling, rotation and translation 
	to a given geometry.
	"""
	# parse the command-line
	args = read_inputs()

	# create the body reading the coordinate file
	body = Geometry(args.coordinates_path)
	
	# re-discretize the body
	body.discretize(n=args.n, cl=args.cl)
	
	# scale the body 
	body.scale(ratio=args.scale)
	# rotation of the body
	body.rotation(aoa=args.rotation, x_rot=args.center[0], y_rot=args.center[1])
	# translation of the body
	body.translation(x_trans=args.translation[0], y_trans=args.translation[1])
	
	# write the new coordinates into a file
	outfile_path = '%s/%s.%s' % (args.save_dir, args.name, args.extension)
	body.write(outfile_path)
	
	# plot the body (initial and new configuration)
	if args.show:
		body.plot()


if __name__ == '__main__':
	main()
