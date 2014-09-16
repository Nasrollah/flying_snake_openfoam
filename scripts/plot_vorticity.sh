#!/bin/sh

# file: $FLYING_SNAKE_OPENFOAM/scripts/run_plot_vorticity.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: run a macro in ParaView to plot the vorticity field


CASE_PATH=$1
FLYING_SNAKE_DIR="/home/mesnardo/flying_snake_openfoam"

echo "\nPlot and Save the vorticity field around the flying snake"
echo "\t--> OpenFOAM case: $CASE_PATH\n"

paraFoam -case $CASE_PATH -touch >> $CASE_PATH/log
pvbatch --use-offscreen-rendering $FLYING_SNAKE_DIR/scripts/plot_vorticity.py \
		--case $CASE_PATH >> $CASE_PATH/log
