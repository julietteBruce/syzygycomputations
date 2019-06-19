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
infile_unzipped=/tmp/$(basename -s .gz $infile)
/bin/zcat -f $infile > $infile_unzipped

## compute rank
export p=$p
export matrixFile=$infile_unzipped

cat $rankScript | /usr/local/bin/gap  -q 

rm $infile_unzipped
exit 
