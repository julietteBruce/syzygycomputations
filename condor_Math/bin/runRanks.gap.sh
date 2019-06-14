#!/bin/bash

infile=$1   # multidegree_4_7_28.dat
rankScript=$2
### if p = 0  RationalField() for the rank calculation;  otherwise GF(p) 
if [ $# -eq 3 ]
then
    p=$3
else
    ## maybe the default should be p = 32003
    p=0
fi

## unzip $infile into /tmp
#infile_unzipped=/tmp/$(basename -s .gz $infile)
#/bin/zcat -f $infile > $infile_unzipped

## compute rank

#/usr/local/bin/magma -b file:=$infile_unzipped p:=$p $rankScript

(export p=32003;  export matrixFile=$infile; cat ../src/ranks.gap | gap -q 1> rank.out 2> rank.err ) &


## this probably doesn't work because MAGMA might not follow exit status conventions...
#magmaStatus=$?

## clean up
#rm $infile_unzipped
#exit $magmaStatus




