#!/bin/sh

# file: $FLYING_SNAKE_OPENFOAM/scripts/run_case.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: run a OpenFOAM simulation given the case path


# path of the case
CASE_PATH=$1

# domain decomposition
decomposePar > simu.log

# run the simulation in parallel
mpirun -np 6 icoFoam -case $CASE_PATH -parallel >> simu.log

# reconstruct the solution
reconstructPar -case $CASE_PATH >> simu.log
