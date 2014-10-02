#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/generate_box_obj.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate a triangulated 2D box and export as .OBJ file


import argparse

import numpy


class Vertex:
	def __init__(self, x ,y, index):
		self.x, self.y = x, y
		self.index = index


class Face:
	def __init__(self, vertex_1, vertex_2, vertex_3):
		self.vertex_1 = vertex_1
		self.vertex_2 = vertex_2
		self.vertex_3 = vertex_3


def main():
	# coorfinates of the 2D box
	x_bl, y_bl = 0.0, 0.0
	x_tr, y_tr = 3.0, 3.0

	# number of points in each directions
	nx, ny = 4, 4

	x = numpy.linspace(x_bl, x_tr, nx)
	y = numpy.linspace(y_bl, y_tr, ny)

	vertices = numpy.empty(nx*ny, dtype=object)
	for j in xrange(ny):
		for i in xrange(nx):
			vertices[j*nx+i] = Vertex(x[i], y[j], j*nx+i)
	
	horizontal_edges = numpy.empty((nx-1)*ny, dtype=object)
	vertical_edges = numpy.empty(nx*(ny-1), dtype=object)
	diagonal_edges = numpy.empty((nx-1)*(ny-1), dtype=object)

	for j in xrange(ny):
		for i in xrange(nx-1):
			horizontal_edges[j*(nx-1)+i] = Edge(vertices[j*nx+i], 
												vertices[j*nx+i+1])

	for j in xrange(ny-1):
		for i in xrange(nx):
			vertical_edges[j*nx+i] = Edge(vertices[j*nx+i], 
										  vertices[(j+1)*nx+i])

	for j in xrange(ny-1):
		for i in xrange(nx-1):
			diagonal_edges[j*(nx-1)+i] = Edge(vertices[j*nx+i+1],
											  vertices[(j+1)*nx+i])


	for i, edge in enumerate(diagonal_edges):
		print '%d - (%g,%g) - (%g,%g)' % (i, 
										  edge.start.x, edge.start.y,
										  edge.end.x, edge.end.y)

	faces = numpy.empty(2*(nx-1)*(ny-1), dtype=object) 

	lower_faces = numpy.empty((nx-1)*(ny-1), dtype=object)
	upper_faces = numpy.empty((nx-1)*(ny-1), dtype=object)

	for j in xrange(ny-1):
		for i in xrange(nx-1):
			lower_faces[j*(nx-1)+i] = Face(horizontal_edges[j*(nx-1)+i],
										   vertical_edges[j*nx+i],
										   diagonal_edges[j*(nx-1)+i])
			upper_faces[j*(nx-1)+i] = Face(horizontal_edges[(j+1)*(nx-1)+i],
										   vertical_edges[j*nx+i+1],
										   diagonal_edges[j*(nx-1)+i])
		

	faces = numpy.insert(upper_faces, 
						 numpy.arange(len(lower_faces)), 
						 lower_faces)


	for i, face in enumerate(upper_faces):
		print 'edge1'
		print '%d - (%g,%g) - (%g,%g)' % (i,
										  face.edge_1.start.x,
										  face.edge_1.start.y,
										  face.edge_1.end.x,
										  face.edge_1.end.y)
		print 'edge2'
		print '%d - (%g,%g) - (%g,%g)' % (i,
										  face.edge_2.start.x,
										  face.edge_2.start.y,
										  face.edge_2.end.x,
										  face.edge_2.end.y)
		print 'edge3'
		print '%d - (%g,%g) - (%g,%g)' % (i,
										  face.edge_3.start.x,
										  face.edge_3.start.y,
										  face.edge_3.end.x,
										  face.edge_3.end.y)

		header = ( '# Wavefront OBJ file\n'
				   '# points: %d\n'
				   '# faces: %d\n'
				   '# zones: 1\n'
				   '# regions: 0 %s\n'
				   % (nx*ny, 2*(nx-1)*(ny-1), 'test'))

		with open('test.obj', 'r') as outfile:
			outfile.write(header)
			for vertex in vertices:
				outfile.write('v %.6f %.6f %.6f\n' % (vertex.x, vertex.y, 1.0))
			outfile.write('g %s\n' % 'test')
			for face in faces:
				outfile.write('f %d %d %d\n' % (face.edge_1.start))


if __name__ == '__main__':
	main()
