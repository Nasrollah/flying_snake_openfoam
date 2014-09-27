#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_pressure_coefficient.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Plot the surface pressure coefficient


import argparse

import numpy
from matplotlib import pyplot


def main():
	infile_path = '/home/mesnardo/tmp/test0.csv'
	with open(infile_path, 'r') as infile:
		p, x = numpy.loadtxt(infile, usecols=(0, 4),dtype=float, 
							 delimiter=',', skiprows=1, unpack=True)
	
	pyplot.figure(figsize=(6, 6))
	pyplot.grid(True)
	pyplot.xlabel(r'$x$', fontsize=18)
	pyplot.ylabel(r'$y$', fontsize=18)
	pyplot.scatter(x, p)
	pyplot.show()


if __name__ == '__main__':
	main()
