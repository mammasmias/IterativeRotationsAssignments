#!/bin/bash


# This script calls the run_single.sh script 50 times for each input file.
# The argument passed is the path to a specific data/.../xyz folder. This 
# script then takes each file in that folder as input to run_single.sh.
#
# Example run command:
#
#   ./run_many.sh data/Al_N/xyz
#
# Before you run, make sure all software is compiled! SEE README


# path passed
p=$1

# for each file in this path
for f in $(ls -v ${p});
do

  # do 50 run_single.sh commands
  for i in {1..50};
  do

    ./run_single.sh $p/$f

  done
done
