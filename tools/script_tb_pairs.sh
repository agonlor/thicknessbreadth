#!/bin/sh

######################################################################
# Script Name: script_tb_pairs
# Description: Take a mesh, convert it into a PGM file, compute its
# 						  thickness-breadth pairs, etc
# Args       : The mesh file name and the size (integer)
# Author     : Aldo Gonzalez-Lorenzo                                                
# Email      : aldo.gonzalez-lorenzo@univ-amu.fr                                          
######################################################################

# -- Arguments
if [ $# -ne 2 ]
then
  echo "Usage: $0 <mesh file> <precision>"
  exit 1
fi
echo "script_tb_pairs.sh: file = "$1".* precision = "$2

# -- Create the workspace
mkdir $1$2
cp ../data/mesh/$1.* $1$2
cd $1$2
curr_dir=$(pwd)

# -- Voxelize the mesh
../binvox -d $2 -v -fit $1.*				# mesh -> binvox
../binvox2pgm $1.binvox $1$2.pgm		# binvox -> pgm
../pgm2obj $1$2.pgm $1$2.obj				# pgm -> obj
rm $1.binvox

# -- Compute the measures
echo "  == Computing the measures =="
../../build/tb $1$2.pgm
../tb-diag.py $1$2.json # make TB diagram
../tb-balls.py $1$2.json 	# make TB balls



# -- Get the other things: balls, small generators, openings and closings
#echo "  == Computing balls =="
#../tools/aldom balls $1$2.pgm
#echo "  == Computing generators =="
#../tools/aldom makereduction $1$2.pgm
#../tools/aldom smallgen $1$2.pgm
#../tools/aldom printcubes $1$2.pgm homgen -r 0.1 -s 0.01
#../tools/aldom printcubes $1$2.pgm cohomgen -r 0.1 -s 0.01
#echo "  == Computing openings and closings =="
#../tools/aldom closeopen $1$2.pgm
#../tools/aldom printcubes $1$2.pgm close -r 0.1 -s 0.01
#../tools/aldom printcubes $1$2.pgm open -r 0.1 -s 0.01
#../tools/aldom printcubes $1$2.pgm complex -r 0.1 -s 0.01