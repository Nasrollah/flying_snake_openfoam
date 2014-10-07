#!/bin/sh

# file: $FLYING_SNAKE_OPENFOAM
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Generate a 2D extruded mesh with GMSH ready for OpenFOAM


# case directory
CASE_PATH=$FLYING_SNAKE_OPENFOAM/test/cylinder_550_xxK_141007_1

# coordinates file name and path
COORD_NAME=cylinder_0.5_10.dat
COORD_PATH=$CASE_PATH/$COORD_NAME

# .geo file name and path
GEO_NAME=cylinder.geo
GEO_PATH=$CASE_PATH/$GEO_NAME

# .msh file name and path
MSH_NAME=cylinder.msh
MSH_PATH=$CASE_PATH/$MSH_NAME

# patch name
BODY_NAME=cylinder

# generate the .geo file
echo generating .geo file...
python generate_geo_file.py --coordinates $COORD_PATH --save-dir $CASE_PATH

# run GMSH
echo running GMSH...
gmsh $GEO_PATH -3 -bgm bgmesh.pos -format msh

# convert mesh to OpenFOAM
echo converting to OpenFOAM...
python gmsh_to_foam.py --case $CASE_PATH --mesh $MSH_PATH --body-name $BODY_NAME

# open ParaFOAM
echo opening ParaFOAM...
paraFoam -case $CASE_PATH
