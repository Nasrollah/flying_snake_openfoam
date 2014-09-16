#!/usr/bin/env python

# file: $FLYING_SNAKE_OPENFOAM/scripts/plot_vorticity.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: macro to run ParaView and plot the vorticity


import argparse
import os

try:
	paraview.simple
except:
	from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the vorticity field '
												 'with ParaFOAM')
	# fill the parser with arguments
	parser.add_argument('--case', dest='case', type=str, default='.',
						help='path of the OpenFOAM case')
	return parser.parse_args()


def main():
	"""Executes the ParaView macro to plot vorticity field."""
	# parse the command-line
	args = read_inputs()
	

	# create images folder if does not exist
	images_path = '%s/images' % args.case
	if not os.path.isdir(images_path):
		os.makedirs(images_path)

	of_file_name = '%s.OpenFOAM' % os.path.basename(os.path.normpath(args.case))

	flying_snake = PV3FoamReader(FileName=('%s/%s' % (args.case, of_file_name)))
	flying_snake.VolumeFields = ['p', 'U']
	flying_snake.MeshParts = ['front - patch']

	# create a render view
	render_view = GetRenderView()
	render_view.CenterAxesVisibility = 0
	render_view.OrientationAxesVisibility = 0
	render_view.CameraPosition = [2.25, 0.0, 8.0]
	render_view.CameraFocalPoint = [2.25, 0.0, 1.0]
	render_view.CameraClippingRange = [6.93, 7.105]
	render_view.CenterOfRotation = [0.0, 0.0, 1.0]
	render_view.CameraParallelScale = 30.0

	# create an animation scene
	animation_scene = GetAnimationScene()

	data_representation = Show()
	data_representation.EdgeColor = [0.0, 0.0, 0.5]	# HSV ?

	compute_derivatives = ComputeDerivatives()
	compute_derivatives.Scalars = ['POINTS', 'p']
	compute_derivatives.Vectors = ['POINTS', 'U']
	compute_derivatives.OutputTensorType = 'Nothing'
	compute_derivatives.OutputVectorType = 'Vorticity'

	data_representation_2 = Show()
	data_representation_2.EdgeColor = [0.0, 0.0, 0.5]	# HSV ?

	text = Text()

	a3_Vorticity_PVLookupTable = GetLookupTableForArray('Vorticity', 3, 
					RGBPoints=[-5.0, 0.0, 0.0, 1.0, 5.0, 1.0, 0.0, 0.0], 
					VectorMode='Component', 
					VectorComponent=2, 
					NanColor=[0.4980392156862745, 0.4980392156862745, 
							  0.4980392156862745],
					ColorSpace='Diverging', 
					ScalarRangeInitialized=1.0, 
					LockScalarRange=1)

	a3_Vorticity_PiecewiseFunction = CreatePiecewiseFunction()

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
	
	GetRenderView().Representations.append(scalar_bar_widget_representation)

	data_representation.Visibility = 0

	data_representation_2.ColorArrayName = 'Vorticity'
	data_representation_2.LookupTable = a3_Vorticity_PVLookupTable
	data_representation_2.ColorAttributeType = 'CELL_DATA'

	text.Text = 'time = 0'

	data_representation_3 = Show()
	data_representation_3.FontSize = 12
	data_representation_3.TextScaleMode = 2
	data_representation_3.Position = [0.02, 0.9]
	data_representation_3.Color = [0.0, 0.0, 0.0]

	view = GetActiveView()
	view.ViewSize = [975, 483]
	
	time_steps = flying_snake.TimestepValues
	for ite in time_steps:
		render_view.ViewTime = ite
		animation_scene.AnimationTime = ite
		text.Text = 'time = %g' % ite
		WriteImage('%s/vorticity_%g.png' % (images_path, ite))


if __name__ == '__main__':
	main()
