#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_geo_file.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate .geo file to be read by GMSH


import argparse

import numpy


class Entity:
	"""Contains info about either the domain of the body."""
	def __init__(self, x, y, cl, counter):
		"""Creates the entity.
		
		Arguments
		---------
		x, y -- coordinates of the entity.
		cl -- characteristic length for the discretization.
		counter -- counters to keep track of number of GMSH entities.
		"""
		self.n = x.size
		self.back = {'points': numpy.empty(self.n, dtype=object),
					 'lines': numpy.empty(self.n, dtype=object),
					 'line-loop': None}
		self.front = {'points': numpy.empty(self.n, dtype=object),
					  'lines': numpy.empty(self.n, dtype=object),
					  'line-loop': None}
		self.back2front = {'lines': numpy.empty(self.n, dtype=object),
						   'line-loops': numpy.empty(self.n, dtype=object)}

		self.create_points(x, y, cl, counter)
		self.create_lines(counter)
		self.create_line_loops(counter)

	def create_points(self, x, y, cl, counter):
		"""Fills the array of points with Point objects.
		
		Arguments
		---------
		x, y -- coordinates of the entity.
		cl -- characteristic length for the discretization.
		counter -- counters to keep track of number of GMSH entities.
		"""
		for i in xrange(self.n):
			self.back['points'][i] = Point(x[i], y[i], 0.0, cl, 
										   counter.points())
		for i in xrange(self.n):
			self.front['points'][i] = Point(x[i], y[i], 1.0, cl,
											counter.points())
	
	def create_lines(self, counter):
		"""Fills the array of lines with Line objects.
		
		Arguments
		---------
		counter -- counters to keep track of number of GMSH entities.
		"""
		for i in xrange(self.n-1):
			self.back['lines'][i] = Line(self.back['points'][i].index, 
										 self.back['points'][i+1].index,
										 counter.lines())
		self.back['lines'][-1] = Line(self.back['points'][-1].index,
									  self.back['points'][0].index,
									  counter.lines())
		for i in xrange(self.n-1):
			self.front['lines'][i] = Line(self.front['points'][i].index, 
										  self.front['points'][i+1].index,
										  counter.lines())
		self.front['lines'][-1] = Line(self.front['points'][-1].index,
									   self.front['points'][0].index,
									   counter.lines())
		for i in xrange(self.n):
			self.back2front['lines'][i] = Line(self.back['points'][i].index,
											   self.front['points'][i].index,
											   counter.lines())

	def create_line_loops(self, counter):
		"""Creates line-loops for different sides.
		
		Arguments
		---------
		counter -- counters to keep track of number of GMSH entities.
		"""
		self.back['line-loop'] = LineLoop([line.index 
										   for line in self.back['lines']],
										  counter.line_loops())
		self.front['line-loop'] = LineLoop([line.index
											for line in self.front['lines']],
										   counter.line_loops())
		for i in xrange(self.n-1):
			self.back2front['line-loops'][i] = LineLoop(
										[self.back['lines'][i].index,
										 self.back2front['lines'][i+1].index,
										 -self.front['lines'][i].index,
										 -self.back2front['lines'][i].index],
										counter.line_loops())
		self.back2front['line-loops'][-1] = LineLoop(
										[self.back['lines'][-1].index,
										 self.back2front['lines'][0].index,
										 -self.front['lines'][-1].index,
										 -self.back2front['lines'][-1].index],
										 counter.line_loops())


class Point:
	"""Contains info about a GMSH point."""
	def __init__(self, x, y ,z, cl, index):
		"""Initializes the point with its coordinates, its characteristic length
		and its identification number.

		Arguments
		---------
		x, y, z -- coordinates of the point.
		cl -- characteristic length of the point.
		index -- identification number of the point.
		"""
		self.index = index
		self.x, self.y, self.z = x, y, z
		self.cl = cl


class Line:
	"""Contains info about a GMSH line."""
	def __init__(self, point_start, point_end, index):
		"""Defines the two ending-points and the identification number.
		
		Arguments
		---------
		point_start, point_end -- identification number of points defining line.
		index -- identification number of the line.
		"""
		self.index = index
		self.start = point_start
		self.end = point_end


class LineLoop:
	"""Contains info about a GMSH line loop."""
	def __init__(self, lines, index):
		"""Stores the identification number of the line loop
		and the lines that compose the loop.
		
		Arguments
		---------
		lines -- idenification number of lines defining the loop.
		index -- identification number of the line loop.
		"""
		self.index = index
		self.lines = lines


class PlaneSurface:
	"""Contains info about a GMSH plane surface."""
	def __init__(self, line_loops, index):
		"""Stores the indentification of the plane surface and the line-loops
		defining the plane surface.
	
		Arguments
		---------
		line_loops -- array of line-loops.
		index -- ideintification number of the plane surface.
		"""
		self.index = index
		self.line_loops = line_loops


class Counter:
	"""Keeps track of number of GMSH entities."""
	def __init__(self):
		"""Initializes counters to zero."""
		self.n_points = 0
		self.n_lines = 0
		self.n_line_loops = 0
		self.n_surfaces = 0
		self.n_surface_loops = 0
		self.n_volumes = 0

	def points(self):
		"""Increments points counter by 1 and returns the value."""
		self.n_points += 1
		return self.n_points

	def lines(self):
		"""Increments lines counter by 1 and returns the value."""
		self.n_lines += 1
		return self.n_lines

	def line_loops(self):
		"""Increments line-loops counter by 1 and returns the value."""
		self.n_line_loops += 1
		return self.n_line_loops

	def surfaces(self):
		"""Increments plane-surfaces counter by 1 and returns the value."""
		self.n_surfaces += 1
		return self.n_surfaces

	def surface_loops(self):
		"""Increments surface-loops counter by 1 and returns the value."""
		self.n_surface_loops += 1
		return self.n_surface_loops

	def volumes(self):
		"""Increments volumes counter by 1 and returns the value."""
		self.n_volumes += 1
		return self.n_volumes


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generate a .geo file',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--coordinates', dest='coordinates_path', type=str,
						help='path of the coordinates file')
	parser.add_argument('--geo-name', dest='geo_name', type=str, 
						default='cylinder',
						help='name of the .geo file (without extension)')
	parser.add_argument('--body-name', dest='body_name', type=str,
						default='cylinder',
						help='name of the body to be used in OpenFOAM')
	parser.add_argument('--save-dir', dest='save_dir', type=str, default='.',
						help='directory where to save .geo file')
	parser.add_argument('--bottom-left', '-bl', dest='bl', type=float,
						nargs='+', default=[-20.0, -20.0],
						help='coordinates of the bottom-left corner of the '
							 'computational domain')
	parser.add_argument('--top-right', '-tr', dest='tr', type=float,
						nargs='+', default=[20.0, 20.0],
						help='coordinates of the top-right corner of the '
							 'computational domain')
	parser.add_argument('--n-exterior', dest='n_exterior', type=int, default=20,
						help='number of points on each external boundaries '
							 '(inlet, outlet, bottom and top)')
	parser.add_argument('--n-segment', dest='n_segment', type=int, default=2,
						help='number of points on each segment of the geometry')
	return parser.parse_args()


def main():
	"""Generates a .geo file that will be read by GMSH to generate the mesh."""
	# parse the command-line
	args = read_inputs()

	# domain coordinates
	x_domain = numpy.array([args.bl[0], args.tr[0], args.tr[0], args.bl[0]])
	y_domain = numpy.array([args.bl[1], args.bl[1], args.tr[1], args.tr[1]])
	bl_x, bl_y = args.bl[0], args.bl[1]	# coordinates of bottom-left corner
	tr_x, tr_y = args.tr[0], args.tr[1]	# coordinates of top-right corner

	# read the coordinates file
	with open(args.coordinates_path, 'r') as infile:
		x, y = numpy.loadtxt(infile, dtype=float, delimiter='\t', skiprows=1, 
							 unpack=True)
	n = x.size	# number of points on the geometry

	print '\ninput\n-----'
	print 'coordinates path: %s' % args.coordinates_path
	print 'number of points: %d' % n

	# compute length of each body-segments
	lengths = numpy.append(numpy.sqrt((x[:-1]-x[1:])**2+(y[:-1]-y[1:])**2),
						   numpy.sqrt((x[0]-x[-1])**2+(y[0]-y[-1])**2))

	print 'minimum segment-length: %g' % lengths.min()
	print 'maximum segment-length: %g' % lengths.max()
	print 'average segment-length: %g' % (lengths.sum()/lengths.size)

	# calculate the characteristic lengths
	
	cl_exterior = (args.tr[0]-args.bl[0])/args.n_exterior
	cl_segment = lengths.max()/args.n_segment

	print '\ncharacteristic lengths\n----------------------'
	print 'external boundaries: %g' % cl_exterior
	print 'geometry (maximum): %g' % cl_segment

	# create and intialize counter to zero
	counter = Counter()

	# create the body
	body = Entity(x, y, cl_segment, counter)

	# create the external domain
	domain = Entity(x_domain, y_domain, cl_exterior, counter)

	surfaces = {'back': [body.back['line-loop'].index, 
						 domain.back['line-loop'].index],
			    'front': [body.front['line-loop'].index,
			   			  domain.front['line-loop'].index]}


	with open('%s/%s.geo' % (args.save_dir, args.geo_name), 'w') as outfile:
		# write characteristic lengths
		outfile.write('// characteristic lengths\n')
		outfile.write('cl__1 = %f;\n' % cl_segment)
		outfile.write('cl__2 = %f;\n' % cl_exterior)
		# write body points
		outfile.write('// body points\n')
		for point in body.back['points']:
			outfile.write('Point(%d) = {%f, %f, %f, %f};\n' 
						  % (point.index, point.x, point.y, point.z, point.cl))
		for point in body.front['points']:
			outfile.write('Point(%d) = {%f, %f, %f, %f};\n' 
						  % (point.index, point.x, point.y, point.z, point.cl))
		# write domain points
		outfile.write('// domain points\n')
		for point in domain.back['points']:
			outfile.write('Point(%d) = {%f, %f, %f, %f};\n' 
						  % (point.index, point.x, point.y, point.z, point.cl))
		for point in domain.front['points']:
			outfile.write('Point(%d) = {%f, %f, %f, %f};\n' 
						  % (point.index, point.x, point.y, point.z, point.cl))
		# write body lines
		outfile.write('// body lines\n')
		for line in body.back['lines']:
			outfile.write('Line(%d) = {%d, %d};\n'
						  % (line.index, line.start, line.end))
		for line in body.front['lines']:
			outfile.write('Line(%d) = {%d, %d};\n'
						  % (line.index, line.start, line.end))
		for line in body.back2front['lines']:
			outfile.write('Line(%d) = {%d, %d};\n'
						  % (line.index, line.start, line.end))
		# write domain lines
		outfile.write('// domain lines\n')
		for line in domain.back['lines']:
			outfile.write('Line(%d) = {%d, %d};\n'
						  % (line.index, line.start, line.end))
		for line in domain.front['lines']:
			outfile.write('Line(%d) = {%d, %d};\n'
						  % (line.index, line.start, line.end))
		for line in domain.back2front['lines']:
			outfile.write('Line(%d) = {%d, %d};\n'
						  % (line.index, line.start, line.end))
		# write body line-loops
		outfile.write('// body line-loops\n')
		outfile.write('Line Loop(%d) = {%s};\n'
					  % (body.back['line-loop'].index, 
					  	 ', '.join([str(i) 
						 			for i in body.back['line-loop'].lines])))
		outfile.write('Line Loop(%d) = {%s};\n'
					  % (body.front['line-loop'].index, 
					  	 ', '.join([str(i) 
						 			for i in body.front['line-loop'].lines])))
		for line_loop in body.back2front['line-loops']:
			outfile.write('Line Loop(%d) = {%s};\n'
					  % (line_loop.index, 
					  	 ', '.join([str(i) for i in line_loop.lines])))
		# write domain line-loops
		outfile.write('// domain line-loops\n')
		outfile.write('Line Loop(%d) = {%s};\n'
					  % (domain.back['line-loop'].index,
					  	 ', '.join([str(i) 
						 			for i in domain.back['line-loop'].lines])))
		outfile.write('Line Loop(%d) = {%s};\n'
					  % (domain.front['line-loop'].index,
					  	 ', '.join([str(i) 
						 			for i in domain.front['line-loop'].lines])))
		for line_loop in domain.back2front['line-loops']:
			outfile.write('Line Loop(%d) = {%s};\n'
					  % (line_loop.index, 
					  	 ', '.join([str(i) for i in line_loop.lines])))
		
		physical_surfaces = {}

		# write plane surfaces
		outfile.write('// plane surfaces\n')
		outfile.write('Plane Surface(%d) = {%d, %d};\n' 
					  % (counter.surfaces(), 
					  	 domain.back['line-loop'].index, 
						 body.back['line-loop'].index))
		physical_surfaces['back'] = counter.n_surfaces
		outfile.write('Plane Surface(%d) = {%d, %d};\n' 
					  % (counter.surfaces(), 
					  	 domain.front['line-loop'].index, 
						 body.front['line-loop'].index))
		physical_surfaces['front'] = counter.n_surfaces
		# write body ruled surfaces
		outfile.write('// body rules surfaces\n')
		physical_surfaces['body'] = []
		for line_loop in body.back2front['line-loops']:
			outfile.write('Ruled Surface(%d) = {%d};\n' 
						  % (counter.surfaces(), line_loop.index))
			physical_surfaces['body'].append(counter.n_surfaces)
		# write domain ruled surfaces
		outfile.write('// domain ruled surfaces\n')
		physical_surfaces['external'] = []
		for line_loop in domain.back2front['line-loops']:
			outfile.write('Ruled Surface(%d) = {%d};\n' 
						  % (counter.surfaces(), line_loop.index))
			physical_surfaces['external'].append(counter.n_surfaces)

		# define surface loop
		outfile.write('// surface-loop\n')
		outfile.write('Surface Loop(%d) = {%s};\n'
					  % (counter.surface_loops(),
					  	 ', '.join([str(i) 
						 			for i in xrange(1, counter.n_surfaces+1)])))
		outfile.write('// volume\n')
		outfile.write('Volume(%d) = {%d};\n' 
					  % (counter.volumes(), counter.n_surface_loops))

		# write physical surfaces
		outfile.write('// physical surfaces\n')
		outfile.write('Physical Surface("back") = {%d};\n' 
					  % physical_surfaces['back'])
		outfile.write('Physical Surface("front") = {%d};\n'
					  % physical_surfaces['front'])
		outfile.write('Physical Surface("inlet") = {%d};\n' 
					  % physical_surfaces['external'][-1])
		outfile.write('Physical Surface("outlet") = {%d};\n'
					  % physical_surfaces['external'][1])
		outfile.write('Physical Surface("bottom") = {%d};\n' 
					  % physical_surfaces['external'][0])
		outfile.write('Physical Surface("top") = {%d};\n'
					  % physical_surfaces['external'][2])
		outfile.write('Physical Surface("%s") = {%s};\n'
					  % (args.body_name, 
						', '.join([str(i) for i in physical_surfaces['body']])))
		outfile.write('Physical Volume("internal") = {%d};\n' 
					  % counter.n_volumes)
		
		# create a field box
		'''outfile.write('// field box\n')
		box_bl_x, box_bl_y = -2.0, -2.0
		box_tr_x, box_tr_y = +2.0, +2.0
		outfile.write('Field[1] = Box;\n')
		outfile.write('Field[1].VIn = %f;\n' % cl_segment)
		outfile.write('Field[1].VOut = %f;\n' % cl_exterior)
		outfile.write('Field[1].XMin = %f;\n' % box_bl_x)
		outfile.write('Field[1].XMax = %f;\n' % box_tr_x)
		outfile.write('Field[1].YMin = %f;\n' % box_bl_y)
		outfile.write('Field[1].YMax = %f;\n' % box_tr_y)
		outfile.write('Background Field = 1;\n')
		
		# parameters for GMSH
		outfile.write('// GMSH parameters\n')
		outfile.write('Recombine Surface{1} = 0;\n')
		outfile.write('Mesh.Algorithm = 8;\n')
		outfile.write('Mesh.Smoothing = 100;\n')
		outfile.write('General.ExpertMode = 1;')'''


if __name__ == '__main__':
	main()
