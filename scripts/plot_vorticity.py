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
	parser.add_argument('--vortlim', dest='vort_lim', type=float, default=5.0,
						help='upper limit of zero-symmetric vorticity range')
	parser.add_argument('--times', dest='times', type=float, nargs='+',
						default=None,
						help='range of times to plot (min, max, increment)')
	parser.add_argument('--start', dest='start', type=float, default=None,
						help='starting-time to plot')
	parser.add_argument('--end', dest='end', type=float, default=None,
						help='ending-time to plot')
	parser.add_argument('--view', dest='view', type=str, default='near_wake',
						help='choose a pre-recorder view')
	return parser.parse_args()


def dict_views():
	"""Returns a dictionary of the different available views."""
	views = {}
	views['near_wake'] = {'CameraPosition': [2.25, 0.0, 8.0],
						  'CameraFocalPoint': [2.25, 0.0, 0.0],
						  'CameraClippingRange': [6.93, 7.105],
						  'ViewSize': [900, 400]}
	views['snake'] = {'CameraPosition': [0.2, 0.0, 4.0],
					  'CameraFocalPoint': [0.2, 0.0, 0.0],
					  'CameraClippingRange': [2.97, 3.045],
					  'ViewSize': [850, 450]}
	views['outlet'] = {'CameraPosition': [12.0, 0.0, 20.0],
					   'CameraFocalPoint': [12.0, 0.0, 0.0],
					   'CameraClippingRange': [18.81, 19.26],
					   'ViewSize': [850, 500]}
	return views


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

	# get pre-recorded views
	views = dict_views()
	# select the appropriate view
	v = views[args.view]
	# set up the view
	view = GetRenderView()
	view.CenterAxesVisibility = 0
	view.OrientationAxesVisibility = 0
	view.CameraPosition = v['CameraPosition']
	view.CameraFocalPoint = v['CameraFocalPoint']
	view.CameraClippingRange = v['CameraClippingRange']
	view.CenterOfRotation = [0.0, 0.0, 1.0]
	view.CameraParallelScale = 30.0
	view.ViewSize = v['ViewSize']
	view.Background = [0.34, 0.34, 0.34]

	# compute the vorticity
	compute_derivatives = ComputeDerivatives()
	compute_derivatives.Scalars = ['POINTS', 'p']
	compute_derivatives.Vectors = ['POINTS', 'U']
	compute_derivatives.OutputTensorType = 'Nothing'
	compute_derivatives.OutputVectorType = 'Vorticity'

	# edit color-map
	a3_Vorticity_PVLookupTable = GetLookupTableForArray('Vorticity', 3, 
					RGBPoints=[-args.vort_lim, 0.0, 0.0, 1.0, 
							   +args.vort_lim, 1.0, 0.0, 0.0], 
					VectorMode='Component', 
					VectorComponent=2, 
					NanColor=[0.0, 0.0, 0.0],
					ColorSpace='Diverging', 
					ScalarRangeInitialized=1.0, 
					LockScalarRange=1)

	# add a scalar bar
	scalar_bar_widget_representation = CreateScalarBar(ComponentTitle='', 
										Title='vorticity', 
										Position2=[0.1, 0.5], 
										Enabled=1, 
										LabelFontSize=12, 
										LabelColor=[0.0, 0.0, 0.0],
										LookupTable=a3_Vorticity_PVLookupTable,
										TitleFontSize=12, 
										TitleColor=[0.0, 0.0, 0.0], 
										Position=[0.02, 0.25])
	view.Representations.append(scalar_bar_widget_representation)

	# show the vorticity field
	data_representation = Show()
	data_representation.ColorArrayName = 'Vorticity'
	data_representation.LookupTable = a3_Vorticity_PVLookupTable
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
		i_end = numpy.where(time_steps > args.end)[0][0]
		time_steps = time_steps[:i_end]
	if args.times:
		start, end, every = args.times[0], args.times[1], args.times[2]
		time_steps = numpy.arange(start, end+every, every)
	
	# time-loop to plot and save the vorticity field
	for time_step in time_steps:
		view.ViewTime = time_step
		text.Text = 'time = %g' % time_step
		WriteImage('%s/vorticity_%s_%g.png' 
				   % (images_path, args.view, time_step))


if __name__ == '__main__':
	main()
