#!/opt/OpenFOAM/ThirdParty-2.2.2/platforms/linux64Gcc/paraview-3.12.0/bin/pvbatch

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_vorticity.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: macro to run ParaView and plot the vorticity


import argparse
import os

import numpy
from paraview.simple import *


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the vorticity field '
												 'with ParaFOAM')
	# fill the parser with arguments
	parser.add_argument('--case', dest='case', type=str, default='.',
						help='path of the OpenFOAM case')
	parser.add_argument('--variable', '-v', dest='variable', type=str,
						help='plot specified variable (vorticity or pressure)')
	parser.add_argument('--vortlim', dest='vort_lim', type=float, default=5.0,
						help='upper limit of zero-symmetric vorticity range')
	parser.add_argument('--times', dest='times', type=float, nargs='+',
						default=None,
						help='range of times to plot (min, max, increment)')
	parser.add_argument('--start', dest='start', type=float, default=None,
						help='starting-time to plot')
	parser.add_argument('--end', dest='end', type=float, default=None,
						help='ending-time to plot')
	parser.add_argument('--view', dest='view', type=str, default=None,
						help='pre-recorded view')
	parser.add_argument('--bl', dest='bottom_left', type=float, nargs='+', 
						default=[-2.0,-2.0],
						help='bottom-left corner of the view')
	parser.add_argument('--tr', dest='top_right', type=float, nargs='+', 
						default=[+2.0,+2.0],
						help='top-right corner of the view')
	parser.add_argument('--width', dest='width', type=float, default=600,
						help='width (in pixel) of the image generated')
	return parser.parse_args()


def get_view(view='snake'):
	"""Returns the limits and the width of the plot.
	
	Arguments
	---------
	view -- name of the pre-recorded view (default 'snake').
	"""
	views = {}
	views['snake'] = {'x_bl': -1.0, 'y_bl': -1.0,
					  'x_tr': 1.5, 'y_tr': 1.0,
					  'width': 600}
	views['near_wake'] = {'x_bl': -1.0, 'y_bl': -2.0,
						  'x_tr': 5.0, 'y_tr': 2.0,
						  'width': 1000}
	views['wake'] = {'x_bl': -2.0, 'y_bl': -3.0,
					  'x_tr': 20.0, 'y_tr': 3.0,
					  'width': 1000}
	return ( views[view]['x_bl'], views[view]['y_bl'],
			 views[view]['x_tr'], views[view]['y_tr'],
			 views[view]['width'] )


def main():
	"""Executes the ParaView macro to plot vorticity field."""
	# parse the command-line
	args = read_inputs()

	# create images folder if does not exist
	images_path = '%s/images' % args.case
	if not os.path.isdir(images_path):
		os.makedirs(images_path)


	# display front patch and read pressure and velocity
	of_file_name = '%s.OpenFOAM' % os.path.basename(os.path.normpath(args.case))
	flying_snake = PV3FoamReader(FileName=('%s/%s' % (args.case, of_file_name)))
	flying_snake.VolumeFields = ['p', 'U']
	flying_snake.MeshParts = ['front - patch']

	# get the limits of the plot and the width of the figure
	if not args.view:
		x_bl, y_bl = args.bottom_left[0], args.bottom_left[1]
		x_tr, y_tr = args.top_right[0], args.top_right[1]
		width = args.width
		args.view = 'custom'
	else:
		x_bl, y_bl, x_tr, y_tr, width = get_view(view=args.view)

	x_center, y_center = 0.5*(x_tr+x_bl), 0.5*(y_tr+y_bl)
	coeff = 20./10.
	coeff = 1.
	h = 0.5*(y_tr-y_bl) + coeff
	height = width*(y_tr-y_bl)/(x_tr-x_bl)

	view = GetRenderView()
	view.ViewSize = [width, height]
	Render()
	view.CenterAxesVisibility = 0
	view.OrientationAxesVisibility = 0
	view.CameraPosition = [x_center, y_center, h]
	view.CameraFocalPoint = [x_center, y_center, 0.0]
	view.CameraViewUp = [0.0, 1.0, 0.0]
	view.CenterOfRotation = [0.0, 0.0, 1.0]
	view.CameraViewAngle = 90.0
	view.Background = [0.34, 0.34, 0.34]
	Render()

	if args.variable == 'vorticity':
		# compute the vorticity
		compute_derivatives = ComputeDerivatives()
		compute_derivatives.Scalars = ['POINTS', 'p']
		compute_derivatives.Vectors = ['POINTS', 'U']
		compute_derivatives.OutputTensorType = 'Nothing'
		compute_derivatives.OutputVectorType = 'Vorticity'
		
		# edit color-map
		PVLookupTable = GetLookupTableForArray('Vorticity', 3, 
									RGBPoints=[-args.vort_lim, 0.0, 0.0, 1.0, 
							   			   	   +args.vort_lim, 1.0, 0.0, 0.0], 
									VectorMode='Component', 
									VectorComponent=2, 
									NanColor=[0.0, 0.0, 0.0],
									ColorSpace='Diverging', 
									ScalarRangeInitialized=1.0, 
									LockScalarRange=1)
	elif args.variable == 'pressure':
		# edit color-map
		PVLookupTable = GetLookupTableForArray('p', 1,
											   RGBPoints=[-1.0, 0.0, 0.0, 1.0,
														  +0.5, 1.0, 0.0, 0.0],
											   VectorMode='Magnitude',
											   NanColor=[0.0, 0.0, 0.0],
											   ColorSpace='HSV',
											   ScalarRangeInitialized=1.0,
											   LockScalarRange=1)

	# add a scalar bar
	scalar_bar = CreateScalarBar(ComponentTitle='', 
								 Title=args.variable, 
								 Position2=[0.1, 0.5], 
								 Enabled=1, 
								 LabelFontSize=12, 
								 LabelColor=[0.0, 0.0, 0.0],
								 LookupTable=PVLookupTable,
								 TitleFontSize=12, 
								 TitleColor=[0.0, 0.0, 0.0], 
								 Position=[0.02, 0.25])
	view.Representations.append(scalar_bar)
	# show the  field
	data_representation = Show()
	if args.variable == 'vorticity':
		array_name = 'Vorticity'
	elif args.variable == 'pressure':
		array_name = 'p'
	data_representation.ColorArrayName = array_name
	data_representation.LookupTable = PVLookupTable
	data_representation.ColorAttributeType = 'CELL_DATA'

	# add text to the view
	text = Text()
	data_representation_3 = Show()
	data_representation_3.FontSize = 12
	data_representation_3.TextScaleMode = 2
	data_representation_3.Position = [0.02, 0.9]	# 0.0, 0.0: bottom-left
	data_representation_3.Color = [0.0, 0.0, 0.0]

	# get the time-steps to plot
	time_steps = numpy.array(flying_snake.TimestepValues)
	if args.start:
		i_start = numpy.where(time_steps >= args.start)[0][0]
		time_steps = time_steps[args.start:]
	if args.end:
		i_end = numpy.where(time_steps >= args.end)[0][0]+1
		time_steps = time_steps[:i_end]
	if args.times:
		start, end, every = args.times[0], args.times[1], args.times[2]
		time_steps = numpy.arange(start, end+every, every)
	
	# time-loop to plot and save the vorticity field
	for time_step in time_steps:
		print 'Time: %g' % time_step
		view.ViewTime = time_step
		text.Text = 'time = %g' % time_step
		WriteImage('%s/%s_%s_%g.png' 
				   % (images_path, args.variable,  args.view, time_step))


if __name__ == '__main__':
	main()
