#!/bin/bash

infile=$1
rankScript=$2

## unzip $infile into /tmp
infile_unzipped=/tmp/$(basename -s .gz $infile)
/bin/zcat -f $infile > $infile_unzipped

## compute rank
/usr/local/bin/magma -b file:=$infile_unzipped $rankScript
## this probably doesn't work because MAGMA might not follow exit status conventions...
magmaStatus=$?

## clean up
rm $infile_unzipped
exit $magmaStatus

