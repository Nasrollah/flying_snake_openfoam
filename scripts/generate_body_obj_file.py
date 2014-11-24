#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_body_obj_file.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: generate an .OBJ file readable by SnappyHexMesh


import os
import argparse

import numpy


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates an .OBJ file '
									 'that will be readable by OpenFOAM '
									 'mesh generator: SnappyHexMesh',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--infile', dest='infile_path', type=str,
						help='path of the coordinates file to be converted')
	parser.add_argument('--name', dest='name', type=str,
						help='name of the .OBJ file generated (no extension)')
	parser.add_argument('--save-dir', dest='save_directory', type=str, 
						default=os.getcwd(),
						help='directory where to save the .obj file')
	return parser.parse_args()


def main():
	"""Generates an .OBJ file from a given coordinates file."""
	# parse the command-line
	args = read_inputs()

	# read the coordinates file
	with open(args.infile_path, 'r') as infile:
		x, y = numpy.loadtxt(infile, dtype=float, 
							 delimiter='\t', skiprows=1, unpack=True)

	# append first element to the end of the array
	tol = 1.0E-06
	if abs(x[0]-x[-1]) > tol and abs(y[0]-y[-1]) > tol:
		x, y = numpy.append(x, x[0]), numpy.append(y, y[0])

	# write .OBJ file
	if not args.name:
		args.name = os.path.basename(os.path.splitext(args.infile_path)[0])
	outfile_path = '%s/%s.obj' % (args.save_directory, args.name)
	header = ( '# Wavefront OBJ file\n'
			   '# points: %d\n'
			   '# faces: %d\n'
			   '# zones: 1\n'
			   '# Regions: 0 %s\n'
			   % (2*x.size, 2*x.size, args.name) )
	with open(outfile_path, 'w') as outfile:
		outfile.write(header)
		for i in xrange(x.size):
			for j in xrange(1, -1, -1):
				outfile.write('v %.6f %.6f %g\n' % (x[i], y[i], j))
		outfile.write('g %s\n' % args.name)
		for i in xrange(1, x.size):
			outfile.write('f %d %d %d\n' % (2*i, 2*i-1, 2*i+1))
			outfile.write('f %d %d %d\n' % (2*i+1, 2*(i+1), 2*i))
		outfile.write('f %d %d %d\n' % (2*x.size, 2*x.size-1, 1))
		outfile.write('f %d %d %d\n' % (1, 2, 2*x.size))


if __name__ == '__main__':
	main()
