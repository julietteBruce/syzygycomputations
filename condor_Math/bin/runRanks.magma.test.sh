#!/bin/bash

infile=$1
rankScript=$2
### if p = 0  RationalField() for the rank calculation;  otherwise GF(p) 
if [ $# -ge 3 ]
then
    p=$3
else
    ## maybe the default should be p = 32003
    p=0
fi
    

## unzip $infile into /tmp
infile_unzipped=/tmp/$(basename -s .gz $infile)
/bin/zcat -f $infile > $infile_unzipped.$$

## compute rank

if [ $# -eq 4 ]
then
    #/usr/local/bin/magma -b file:=$infile_unzipped.$$ p:=$p transpose:=$4 $rankScript
    /usr/local/bin/magma  file:=$infile_unzipped.$$ p:=$p transpose:=$4 $rankScript
else
    #/usr/local/bin/magma -b file:=$infile_unzipped.$$ p:=$p $rankScript
    /usr/local/bin/magma  file:=$infile_unzipped.$$ p:=$p $rankScript
fi

## this probably doesn't work because MAGMA might not follow exit status conventions...
magmaStatus=$?

## clean up
rm $infile_unzipped.$$
exit $magmaStatus

