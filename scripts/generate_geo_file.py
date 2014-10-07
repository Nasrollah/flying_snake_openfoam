#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_geo_file.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate .geo file to be read by GMSH


import argparse

import numpy


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
						help='name of the patch used in OpenFOAM')
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

	# write .geo file
	with open('%s/%s.geo' % (args.save_dir, args.geo_name), 'w') as outfile:
		# write body points
		outfile.write('// body points\n')
		for i in xrange(n):
			outfile.write('Point(%d) = {%f, %f, 0.0, %f};\n' 
						  % (i+1, x[i], y[i], cl_segment))
		# write body lines
		outfile.write('// body lines\n')
		for i in xrange(n-1):
			outfile.write('Line(%d) = {%d, %d};\n' % (i+1, i+1, i+2))
		outfile.write('Line(%d) = {%d, %d};\n' % (n, n, 1))
		# write domain points
		outfile.write('// domain points\n')
		outfile.write('Point(%d) = {%f, %f, 0.0, %f};\n' 
					  % (n+1, args.bl[0], args.bl[1], cl_exterior))
		outfile.write('Point(%d) = {%f, %f, 0.0, %f};\n' 
					  % (n+2, args.tr[0], args.bl[1], cl_exterior))
		outfile.write('Point(%d) = {%f, %f, 0.0, %f};\n' 
					  % (n+3, args.tr[0], args.tr[1], cl_exterior))
		outfile.write('Point(%d) = {%f, %f, 0.0, %f};\n' 
					  % (n+4, args.bl[0], args.tr[1], cl_exterior))
		# write domain lines
		outfile.write('// domain lines\n')
		outfile.write('Line(%d) = {%d, %d};\n' % (n+1, n+1, n+2))
		outfile.write('Line(%d) = {%d, %d};\n' % (n+2, n+2, n+3))
		outfile.write('Line(%d) = {%d, %d};\n' % (n+3, n+3, n+4))
		outfile.write('Line(%d) = {%d, %d};\n' % (n+4, n+4, n+1))
		# write body line-loop
		outfile.write('// body line-loop\n')
		outfile.write('Line Loop(1) = {%s};\n' 
					  % ', '.join(['%s' % str(i+1) for i in xrange(n)]))
		# write domain line-loop
		outfile.write('// domain line-loop\n')
		outfile.write('Line Loop(2) = {%s};\n' 
					  % ', '.join(['%s' % str(n+i) for i in [1, 2, 3, 4]]))
		# write plane surface
		outfile.write('// plane surface\n')
		outfile.write('Plane Surface(1) = {1, 2};\n')
		# physical surfaces and volume
		outfile.write('// physical surfaces and volume\n')
		outfile.write('Physical Surface("back") = {%d};\n' % 1)
		outfile.write('Physical Surface("front") = {%d};\n' % 86)
		outfile.write('Physical Surface("inlet") = {%d};\n' % 73)
		outfile.write('Physical Surface("outlet") = {%d};\n' % 81)
		outfile.write('Physical Surface("bottom") = {%d};\n' % 85)
		outfile.write('Physical Surface("top") = {%d};\n' % 77)
		outfile.write('Physical Surface("%s") = {%s};\n' 
					  % (args.body_name, 
					  	 ', '.join(str(i) for i in [33, 37, 41, 45, 
						 							49, 53, 57, 61, 65, 69])))
		outfile.write('Physical Volume(1) = {1};\n')
		# create a field box
		outfile.write('// field box\n')
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
		# recombine and extrude to get a 3D mesh with 1 cell in 3rd-direction
		outfile.write('// GMSH parameters\n')
		outfile.write('Recombine Surface{1} = 0;\n')
		outfile.write('Mesh.Algorithm = 8;\n')
		outfile.write('Extrude {0, 0, 1} {'
					  '\nSurface{1};\nLayers{1};\nRecombine;\n'
					  '}\n')

		outfile.write('Mesh.Smoothing = 100;\n')
		outfile.write('General.ExpertMode = 1;\n')


if __name__ == '__main__':
	main()
