#!/bin/bash

#note this is a wrapper arround the matlab script that helps with input/output setup
infile=$1
outdir=${2%/}

outfilename=$(basename ${infile} .dat)

echo "./run_final_qr.sh /usr/local/MATLAB/CURRENT/ ${infile} ${outdir}/${outfilename}.txt"

exec ./run_final_qr.sh /usr/local/MATLAB/CURRENT/ ${infile} ${outdir}/${outfilename}.txt

