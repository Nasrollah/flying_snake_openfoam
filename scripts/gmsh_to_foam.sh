#!/bin/sh

# file: $FLYING_SNAKE_OPENFOAM/scripts/gmsh_to_foam.sh
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: convert a gmsh mesh to an OpenFOAM one


# path of OpenFOAM case
CASE_PATH=$1

# copy the GMSH file inside the case folder
cp $CASE_PATH/mesh_generation/*.msh $CASE_PATH/.

# convert the mesh to a OpenFOAM compatible one
gmshToFoam -case $CASE_PATH *.msh

# change boundary patches in boundary file
python gmsh_to_foam.py --case $CASE_PATH

# check the quality of the mesh
checkMesh -case $CASE_PATH
