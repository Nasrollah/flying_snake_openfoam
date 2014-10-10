#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_box_obj_file.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate a triangulated 2D box and export as .OBJ file


import argparse
import os

import numpy


class Vertex:
	"""Contains info about a vertex."""
	def __init__(self, x ,y, index):
		"""Initializes a vertex by its coordinates and its index in the mesh.
		
		Arguments
		---------
		x, y -- coordinates of the vertex.
		index -- index of the vertex in the mesh.
		"""
		self.x, self.y = x, y
		self.index = index


class Face:
	"""Contains info about a face."""
	def __init__(self, vertex_1, vertex_2, vertex_3):
		"""Initializes a face by its vertices.
		
		Arguments
		---------
		vertex_1, vertex_2, vertex_3 -- Vertex objects that compose the face.
		"""
		self.vertex_1 = vertex_1
		self.vertex_2 = vertex_2
		self.vertex_3 = vertex_3


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates a 2D triangulated '
												 'box and writes into a '
												 'Wavefront OBJ file',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--bottom-left', '-bl', dest='bl', type=float, 
						nargs='+', default=[-2.0, -2.0],
						help='coordinates of the bottom-left corner of the box')
	parser.add_argument('--top-right', '-tr', dest='tr', type=float, 
						nargs='+', default=[2.0, 2.0],
						help='coordinates of the top-right corner of the box')
	parser.add_argument('-z', dest='z', type=float, default=1.0,
						help='z-coordinate of the box')
	parser.add_argument('-n', dest='n', type=int, nargs='+', default=[100, 100],
						help='number of points in the x- and y- directions')
	parser.add_argument('--name', dest='name', type=str, default='box',
						help='name of the OBJ file without the extension')
	parser.add_argument('--save-dir', dest='save_directory', type=str, 
						default=os.getcwd(),
						help='directory where to save the .obj file')
	return parser.parse_args()


def main():
	"""Generates a 2D triangulated box and writes in Wavefront OBJ file."""
	# parse the command-line
	args = read_inputs()

	# coorfinates of the 2D box
	x_bl, y_bl = args.bl[0], args.bl[1]
	x_tr, y_tr = args.tr[0], args.tr[1]

	# number of points in each directions
	nx, ny = args.n[0], args.n[1]

	# create x-stations and y-stations
	x = numpy.linspace(x_bl, x_tr, nx)
	y = numpy.linspace(y_bl, y_tr, ny)

	# store the vertices
	vertices = numpy.empty(nx*ny, dtype=object)
	for j in xrange(ny):
		for i in xrange(nx):
			vertices[j*nx+i] = Vertex(x[i], y[j], j*nx+i)


	# create the faces (lower and upper triangles)
	lower_faces = numpy.empty((nx-1)*(ny-1), dtype=object)
	upper_faces = numpy.empty((nx-1)*(ny-1), dtype=object)
	for j in xrange(ny-1):
		for i in xrange(nx-1):
			lower_faces[j*(nx-1)+i] = Face(vertices[j*nx+i],
										   vertices[j*nx+i+1],
										   vertices[(j+1)*nx+i])
			upper_faces[j*(nx-1)+i] = Face(vertices[(j+1)*nx+i+1],
										   vertices[(j+1)*nx+i],
										   vertices[j*nx+i+1])
	faces = numpy.insert(upper_faces, 
						 numpy.arange(len(lower_faces)), 
						 lower_faces)

	# write the triangulated box into an OBJ file
	header = ( '# Wavefront OBJ file\n'
			   '# points: %d\n'
			   '# faces: %d\n'
			   '# zones: 1\n'
			   '# regions: 0 %s\n'
			   % (nx*ny, 2*(nx-1)*(ny-1), args.name))
	obj_path = '%s/%s.obj' % (args.save_directory, args.name)
	with open(obj_path, 'w') as outfile:
		outfile.write(header)
		for vertex in vertices:
			outfile.write('v %.6f %.6f %.6f\n' % (vertex.x, vertex.y, args.z))
		outfile.write('g %s\n' % args.name)
		for face in faces:
			outfile.write('f %d %d %d\n' % (face.vertex_1.index+1,
											face.vertex_2.index+1,
											face.vertex_3.index+1))
	

if __name__ == '__main__':
	main()
