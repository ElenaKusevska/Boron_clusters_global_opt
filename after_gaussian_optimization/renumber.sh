#!/bin/bash

# So that the numbers are starting from 1:
j=0
for i in {1..50000}
do
   if [ -f $i.xyz ]; then
      let "j=$j+1"
      echo $i
      echo $j
      cp $i.xyz $j.xyz
      rm $i.xyz
   fi
done
