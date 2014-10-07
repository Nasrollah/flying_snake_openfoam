#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/gmsh_to_foam.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: change the boundary patches in folder constant/polyMesh


import os
import argparse


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Change the boundary patches '
									 'in the file constant/polyMesh/boundary '
									 'of the OpenFOAM case')
	# fill the parser with arguments
	parser.add_argument('--case', dest='case', type=str, default='.',
						help='path of the OpenFOAM case')
	parser.add_argument('--mesh', dest='mesh', type=str, default=None,
						help='path of the GMSH file')
	parser.add_argument('--body-name', dest='body_name', type=str, 
						default='cylinder',
						help='name of the body patch in OpenFOAM')
	return parser.parse_args()


def main():
	"""Changes the boundary patches in the file constant/polyMesh/boundary
	of the OpenFOAM case.
	"""
	# parse the command-line
	args = read_inputs()

	# copy the GMSH file inside the case folder
	if args.mesh == None:
		args.mesh = '%s/mesh_generation/*.msh' % args.case
	os.system('cp %s %s/.' % (args.mesh, args.case))

	mesh_path = '%s/*.msh' % args.case
	log_path = '%s/mesh.log' % args.case

	# run OpenFOAM utility gmshToFoam
	os.system('gmshToFoam -case %s %s > %s' % (args.case, mesh_path, log_path))

	# path of the file containing the patches
	boundary_path = '%s/constant/polyMesh/boundary' % args.case

	# read the boundary file
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
	os.system('checkMesh -case %s >> %s' % (args.case, log_path))

	# delete mesh file
	os.system('rm -rf %s' % mesh_path)


if __name__ == '__main__':
	main()
