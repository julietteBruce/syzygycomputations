#!/bin/bash
~/Downloads/sage-7.4/sage ../src/ranks.sage matricies/map_$1_$2/ | sed 's/[,()]//g' | python3 ../src/condor/ranks_to_betti.py $1 $2
