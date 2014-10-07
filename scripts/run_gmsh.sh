#!/bin/sh

# file: $FLYING_SNAKE_OPENFOAM
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate a 2D extruded mesh with GMSH


# path of the .geo file
GEO_PATH=$1

# run GMSH
gmsh $GEO_PATH -3 -bgm bgmesh.pos -format msh
