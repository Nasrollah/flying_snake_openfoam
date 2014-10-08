#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_geo_file.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate .geo file to be read by GMSH


import argparse
import os

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
	parser.add_argument('--n-exterior', dest='n_exterior', type=int,
						default=20,
						help='number of points on inlet boundary '
							 '(inlet, outlet, bottom and top)')
	parser.add_argument('--cl-exterior', dest='cl_exterior', type=float, 
						help='characteristic-length of elements at external '
							 'boundaries')
	parser.add_argument('--n-segment', dest='n_segment', type=int, default=2,
						help='number of points on each segment of the '
							 'dicretized geometry')
	parser.add_argument('--cl-segment', dest='cl_segment', type=float, 
						help='characteristic-length of elements on the body')
	parser.add_argument('--box', dest='box', type=float, nargs='+',
						help='adds volumes of refinement '
							 '(x_bl, y_bl, x_tr, y_tr, cl_in, cl_max)')
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
	print 'number of body points: %d' % n

	# compute length of each body-segments
	lengths = numpy.append(numpy.sqrt((x[:-1]-x[1:])**2+(y[:-1]-y[1:])**2),
						   numpy.sqrt((x[0]-x[-1])**2+(y[0]-y[-1])**2))

	print 'minimum segment-length: %g' % lengths.min()
	print 'maximum segment-length: %g' % lengths.max()
	print 'average segment-length: %g' % (lengths.sum()/lengths.size)

	# domain parameters
	x_bl, y_bl = args.bl[0], args.bl[1]
	x_tr, y_tr = args.tr[0], args.tr[1]

	# external characteristic-length
	if not args.cl_exterior:
		cl_exterior = (y_tr-y_bl)/args.n_exterior
	else:
		cl_exterior = args.cl_exterior

	# body characteristic-length
	if not args.cl_segment:
		cl_segment = lengths.max()/args.n_segment
	else:
		cl_segment = args.cl_segment

	print '\ncharacteristic lengths\n----------------------'
	print 'external boundaries: %g' % cl_exterior
	print 'geometry (maximum): %g' % cl_segment

	# write .geo file
	geo_path = '%s/%s.geo' % (args.save_dir, args.geo_name)
	with open(geo_path, 'w') as outfile:
		# write characteristic-lengths
		outfile.write('cl_exterior = %f;\n' % cl_exterior)
		outfile.write('cl_segment = %f;\n' % cl_segment)
		# write body points
		outfile.write('// body points\n')
		for i in xrange(n):
			outfile.write('Point(%d) = {%f, %f, 0.0, %f};\n' 
						  % (i+1, x[i], y[i], 2*cl_segment))
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
		# physical volume
		outfile.write('// physical volume\n')
		outfile.write('Physical Volume(1) = {1};\n')
		# create field boxes
		if args.box:
			n_level = 10	# number of cells between levels
			x_bl, y_bl = args.box[0], args.box[1]	# bottom-left corner of box
			x_tr, y_tr = args.box[2], args.box[3]	# top-right corner of box
			# parameters of the box
			box = {'id': 1,
				   'x_bl': args.box[0], 'y_bl': args.box[1],
				   'x_tr': args.box[2], 'y_tr': args.box[3],
				   'cl_in': cl_segment, 'cl_out': cl_exterior}
			# write different boxes
			while box['cl_in'] < box['cl_out']:
				outfile.write('// box field %d\n' % box['id'])
				outfile.write('Field[%d] = Box;\n' % box['id'])
				outfile.write('Field[%d].VIn = %f;\n' 
							  % (box['id'], box['cl_in']))
				outfile.write('Field[%d].VOut = %f;\n' 
							  % (box['id'], box['cl_out']))
				outfile.write('Field[%d].XMin = %f;\n' 
							  % (box['id'], box['x_bl']))
				outfile.write('Field[%d].XMax = %f;\n' 
							  % (box['id'], box['x_tr']))
				outfile.write('Field[%d].YMin = %f;\n' 
							  % (box['id'], box['y_bl']))
				outfile.write('Field[%d].YMax = %f;\n' 
							  % (box['id'], box['y_tr']))
				# parameters of the box at the next level
				box['id'] += 1
				box['cl_in'] *= 2.0
				box['x_bl'] -= (n_level-1)*box['cl_in']
				box['y_bl'] -= (n_level-1)*box['cl_in']
				box['x_tr'] += (n_level-1)*box['cl_in']
				box['y_tr'] += (n_level-1)*box['cl_in']
			# write background field
			outfile.write('// background field\n')
			outfile.write('Field[%d] = Min;\n' % box['id'])
			outfile.write('Field[%d].FieldsList = {%s};\n' 
						  % (box['id'],
						     ', '.join([str(i) for i in range(1,box['id'])])))
			outfile.write('Background Field = %d;\n' % box['id'])
		# recombine and extrude to get a 3D mesh with 1 cell in 3rd-direction
		outfile.write('// GMSH parameters\n')
		outfile.write('Recombine Surface{1} = 0;\n')
		outfile.write('Mesh.Algorithm = 8;\n')	# DelQuad algorithm
		outfile.write('Extrude {0, 0, 1} {'
					  '\nSurface{1};\nLayers{1};\nRecombine;\n'
					  '}\n')
		outfile.write('Mesh.Smoothing = 100;\n')
		outfile.write('General.ExpertMode = 1;\n')


if __name__ == '__main__':
	main()
