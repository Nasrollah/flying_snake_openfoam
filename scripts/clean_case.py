#!/usr/bin/env/ python

# file $FLYING_SNAKE_OPENFOAM/scripts/clean_case.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Clean an OpenFOAM simulation folder


import argparse
import os


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Cleans an OpenFOAM '
												 'simulation folder',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--case', dest='case_directory', type=str, default='.',
						help='directory of the OpenFOAM case')
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
	parser.add_argument('--no-post-processing', dest='post_processing', 
						action='store_false',
						help='does not remove the post-processing folder')
	parser.set_defaults(images=True, processors=True, solutions=True, logs=True,
						post_processing=True)
	return parser.parse_args()


def main():
	"""Cleans an OpenFOAM simulation case."""
	# parse the command-line
	args = read_inputs()

	# store different paths into a dictionary if no flag
	parts = {}
	if args.images:
		parts['images'] = '%s/images' % args.case_directory
	if args.processors:
		parts['processors'] = '%s/processor*' % args.case_directory
	if args.solutions:
		parts['solutions'] = ('%s/[1-9]* %s/0.*' 
							  % (args.case_directory, args.case_directory))
	if args.logs:
		parts['logs'] = '%s/*log*' % args.case_directory
	if args.post_processing:
		parts['post-processing'] = '%s/postProcessing' % args.case_directory

	# remove paths that are in the dictionary
	print 'case directory: %s' % args.case_directory
	for key, part in parts.iteritems():
		print '\t--> deleting %s ...' % key
		os.system('rm -rf %s' % part)


if __name__ == '__main__':
	main()
