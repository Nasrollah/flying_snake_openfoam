#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/run_gmsh.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Run GMSH to generate 2D extruded mesh


import argparse
import os


def read_inputs():
	"""Parses the command-line."""
	# create parser
	parser = argparse.ArgumentParser(description='Runs GMSH to generate '
												 'a 2D extruded mesh',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill parser with arguments
	parser.add_argument('--geo', dest='geo_path', type=str, 
						default='./cylinder.geo',
						help='path of the .geo file')
	return parser.parse_args()


def main():
	"""Runs GMSH to generate a 2D extruded mesh."""
	# parse the command-line
	args = read_inputs()

	# run GMSH in batch mode (-3: 3D generation, -format msh: output format)
	os.system('gmsh %s -3 -format msh' % args.geo_path)


if __name__ == '__main__':
	main()
