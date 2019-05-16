#!/bin/bash


bettiFile=$1 # input file with 
outDir=$2

#python3 src/setup_matrices_prod_Pn.py $outDir $bettiFile

for map in $(find $outDir/matrices/* -type d);
do
    python3 src/compute_rank_magma.py $map $outDir/ranks;
done

for rk in $(find $outDir/ranks/*);
do
    magma< $rk;
done
