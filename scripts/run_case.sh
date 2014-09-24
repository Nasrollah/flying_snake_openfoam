#!/bin/sh

# file: $FLYING_SNAKE_OPENFOAM/scripts/run_case.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: run a OpenFOAM simulation given the case path


# path of the case
CASE_PATH=$1

LOG_PATH=$CASE_PATH/simu.log

# domain decomposition
decomposePar -case $CASE_PATH > $LOG_PATH

# run the simulation in parallel
$OPENFOAM_MPIRUN -np 6 icoFoam -case $CASE_PATH -parallel >> $LOG_PATH

# reconstruct the solution
reconstructPar -case $CASE_PATH -newTimes >> $LOG_PATH
