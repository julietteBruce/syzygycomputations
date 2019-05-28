import numpy as np
import sys
import re
import itertools
from os import listdir
from os.path import isfile, join
import os
import argparse
import glob

argparser = argparse.ArgumentParser();
argparser.add_argument('ranks_dir')
argparser.add_argument('betti_dir')

args = argparser.parse_args()

ranks_dir=args.ranks_dir
betti_dir=args.betti_dir

def filename_to_indices(name):
    """Take a filename of the format *ranks/ranks_* and reads off the multidegree"""
    match = re.search('ranks/ranks((?:_\d+)*)',name)
    return tuple(map(int,re.findall("\d+",match.group(1))))

onlyfiles = glob.glob(os.path.join(ranks_dir,"*.txt"))
pqs = [filename_to_indices(f) for f in onlyfiles]

rks_cols = {}
for pq in pqs:
    rks_cols[pq]={}
    with open(os.path.join(ranks_dir,'ranks_{}_{}.txt'.format(pq[0],pq[1]))) as f:
        for line in f:
            (md, rk_col) = line.split(')')
            rks_cols[pq][eval(md+')')] = tuple(map(int, rk_col[1:].split(' ')))

betti={}
for pq in pqs:
    betti[pq]={}
    outBetti = open(os.path.join(betti_dir,'betti_{}_{}.txt'.format(pq[0],pq[1])),"w+")
    pq1=(pq[0]+1, pq[1]-1)
    for md in rks_cols[pq].keys():
        if pq1[1]<0:
            betti[pq][md] = rks_cols[pq][md][1] - rks_cols[pq][md][0]
        else:
            if md not in rks_cols[pq1].keys():
                betti[pq][md] =  rks_cols[pq][md][1] - rks_cols[pq][md][0]
            else:
                betti[pq][md] = rks_cols[pq][md][1] - rks_cols[pq][md][0] - rks_cols[pq1][md][0]
        outBetti.write('{} {}\n'.format(md,betti[pq][md]))
    outBetti.close()
    outSeries=open(os.path.join(betti_dir,'bettiSeries_{}_{}.txt'.format(pq[0],pq[1])),"w+")
    f="+".join(["{}*t_0^({})*t_1^({})*s_0^({})*s_1^({})".format(betti[pq][md],md[0],md[1],md[2],md[3]) for md in betti[pq].keys()])
    outSeries.write(f)
    outSeries.close()


print(betti[(13,0)])

print(sum([betti[(10,2)][md] for md in betti[(10,2)].keys()] ))  
print(sum([betti[(11,2)][md] for md in betti[(11,2)].keys()] ))    
print(sum([betti[(11,1)][md] for md in betti[(11,1)].keys()] ))
print(sum([betti[(12,1)][md] for md in betti[(12,1)].keys()] ))
