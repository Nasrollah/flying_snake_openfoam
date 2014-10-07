#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/gmsh_to_foam.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Change the boundary patches in folder constant/polyMesh


import os
import argparse


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Change the boundary patches '
									 'in the file constant/polyMesh/boundary '
									 'of the OpenFOAM case',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--case', dest='case_directory', type=str, default='.',
						help='directory of the OpenFOAM case')
	parser.add_argument('--mesh', dest='mesh_path', type=str, default='./*.msh',
						help='path of the GMSH .msh file')
	parser.add_argument('--body-name', dest='body_name', type=str, 
						default='cylinder',
						help='name of the body patch used in OpenFOAM')
	return parser.parse_args()


def main():
	"""Changes the boundary patches in the file constant/polyMesh/boundary
	of the OpenFOAM case.
	"""
	# parse the command-line
	args = read_inputs()

	# log file
	log_path = '%s/mesh.log' % args.case_directory

	# run OpenFOAM utility gmshToFoam
	os.system('gmshToFoam -case %s %s > %s' 
			  % (args.case_directory, args.mesh_path, log_path))

	# path of the file containing the patch names
	boundary_path = '%s/constant/polyMesh/boundary' % args.case_directory

	# read the current boundary file
	with open(boundary_path, 'r') as infile:
		lines = infile.readlines()

	# change boundary patches
	for i, line in enumerate(lines):
		if 'front' in line or 'back' in line:
			lines[i+2] = lines[i+2].replace('patch', 'empty')
		elif 'top' in line or 'bottom' in line:
			lines[i+2] = lines[i+2].replace('patch', 'symmetryPlane')
		elif 'inlet' in line:
			lines[i+3] = lines[i+3].replace('patch', 'inlet')
		elif 'outlet' in line:
			lines[i+3] = lines[i+3].replace('patch', 'outlet')
		elif args.body_name in line:
			lines[i+2] = lines[i+2].replace('patch', 'wall')
			lines[i+3] = lines[i+3].replace('patch', 'wall')
			
	# write the new boundary file
	with open(boundary_path, 'w') as outfile:
		outfile.write(''.join(lines))

	# check the quality of the mesh
	os.system('checkMesh -case %s >> %s' % (args.case_directory, log_path))


if __name__ == '__main__':
	main()
