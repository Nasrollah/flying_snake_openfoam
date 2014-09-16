#!/usr/bin/env/ python

# file $FLYING_SNAKE_OPENFOAM/scripts/clean_case.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: script to clean an OpenFOAM simulation folder


import argparse
import os


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Cleans an OpenFOAM '
												 'simulation folder')
	# fill the parser with arguments
	parser.add_argument('--case', dest='case', type=str, default='.',
						help='path of the OpenFOAM case')
	parser.add_argument('--no-images', dest='images', action='store_false',
						help='does not remove the images folder')
	parser.add_argument('--no-processors', dest='processors', 
						action='store_false',
						help='does not remove processor folders')
	parser.add_argument('--no-solutions', dest='solutions', 
						action='store_false',
						help='does not remove solution folders')
	parser.add_argument('--no-logs', dest='logs', action='store_false',
						help='does not remove log files')
	parser.add_argument('--no-of_files', dest='of_files', action='store_false',
						help='does not remove .OpenFOAM files')
	parser.set_defaults(images=True, processors=True, solutions=True, logs=True,
						of_files=True)
	return parser.parse_args()


def main():
	"""Cleans an OpenFOAM simulation case."""
	# parse the command-line
	args = read_inputs()

	# store different paths into a dictionary if no flag
	parts = {}
	if args.images:
		parts['images'] = '%s/images/' % args.case
	if args.processors:
		parts['processors'] = '%s/processor*/' % args.case
	if args.of_files:
		parts['of_files'] = '%s/*.OpenFOAM' % args.case
	if args.solutions:
		parts['solutions'] = '%s/[1-9]*' % args.case
	if args.logs:
		parts['logs'] = '%s/*log*' % args.case

	print 'case: %s' % args.case
	# remove paths that are in the dictionary
	for key, part in parts.iteritems():
		print '\t--> deleting %s ...' % key
		os.system('rm -rf %s' % part)


if __name__ == '__main__':
	main()
