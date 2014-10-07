#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/write_physical_surfaces.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# desription: Write physical surfaces into .geo file


import argparse

import numpy


def read_inputs():
	"""Parses the command-line."""
	# create parser
	parser = argparse.ArgumentParser(description='Write physical surfaces '
												 'into .geo file',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill parser with arguments
	parser.add_argument('--geo', dest='geo_path', type=str,
						help='path of the .geo file')
	parser.add_argument('--back', dest='back', type=int,
						help='identification number of the back surface')
	parser.add_argument('--front', dest='front', type=int,
						help='identification number of the front surface')
	parser.add_argument('--inlet', dest='inlet', type=int,
						help='identification number of the inlet surface')
	parser.add_argument('--outlet', dest='outlet', type=int,
						help='identification number of the outlet surface')
	parser.add_argument('--bottom', dest='bottom', type=int,
						help='identification number of the bottom surface')
	parser.add_argument('--top', dest='top', type=int,
						help='identification number of the top surface')
	parser.add_argument('--body-name', dest='body_name', type=str, 
						default='cylinder',
						help='name of the OpenFOAM patch for the body')
	parser.add_argument('--body', dest='body', type=int, nargs='+',
						help='starting and ending identification numbers '
							 'followed by the increment')
	return parser.parse_args()


def main():
	"""Writes physical surfaces into the .geo file."""
	# parse the command-line
	args = read_inputs()

	# write physical surfaces into .geo file
	with open('%s' % args.geo_path, 'a') as outfile:
		outfile.write('\n// physical surfaces\n')
		outfile.write('Physical Surface("back") = {%d};\n' % args.back)
		outfile.write('Physical Surface("front") = {%d};\n' % args.front)
		outfile.write('Physical Surface("inlet") = {%d};\n' % args.inlet)
		outfile.write('Physical Surface("outlet") = {%d};\n' % args.outlet)
		outfile.write('Physical Surface("bottom") = {%d};\n' % args.bottom)
		outfile.write('Physical Surface("top") = {%d};\n' % args.top)
		outfile.write('Physical Surface("%s") = {%s};\n' 
					  % (args.body_name, 
					  	 ', '.join([str(i) 
						 			for i in numpy.arange(args.body[0],
														  args.body[1]+1,
														  args.body[2])])))


if __name__ == '__main__':
	main()
