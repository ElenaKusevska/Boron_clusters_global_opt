#!/bin/bash

# This is a script, for
# when the outputs come out of the computational center, to separate
# them in three different folders based on the method that was used
# to prepare the initial structure. Because I wanted to look at each of
# the three methods separately, and see which one gave the best initial
# structures.

# one directory for each method:
# (the names are weird because of fortran format - they all
# had the same number of characters)
mkdir combine_from_grid random_coordinate coordinate_walkin 

# the inputs/outputs for every calculation are in folders named 1, 2, 
# 3, ...
for i in {1..324}
do
   cd $i
   if [ -f $i.xyz ]; then
      j=$(cat $i.xyz | head -2 | tail -1| cut -c 1-17) # where the name of the method was
      echo i $i
      echo j $j
      cd ../
      cp -r $i $j
      trash $i
   else 
      cd ../
   fi
done
   

