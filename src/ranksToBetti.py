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
indices = [filename_to_indices(f) for f in onlyfiles]

for ind in indices:
    d = {}
    with open(os.path.join(ranks_dir,'/ranks_{}_{}.txt'.format(ind[0],ind[1]))) as f:
        for line in f:
            (key, val) = line.split(')')
            d[eval(key+')')] = map(int, val[1:-1].split())
    if ind[1] != 0:
        E = {}
        with open(os.path.join(ranks_dir,'/ranks_{}_{}.txt'.format(ind[0]+1,ind[1]-1))) as f:
            for line in f:
                (key, val) = line.split(')')
                E[eval(key+')')] = map(int, val[1:-1].split())
        betti = {}
        for key in d.keys():
            if key in E:
                betti[key] = ((d[key][1] - d[key][0]) - E[key][0])
            else:
                betti[key] = (d[key][1] - d[key][0])
    else:
        betti = {}
        for key in d.keys():
            betti[key] = (d[key][1] - d[key][0])
    outBetti = open(betti_dir+'/betti_{}_{}.txt'.format(ind[0],ind[1]),"w+")
    outBetti.writelines(['{} {}\n'.format(key,betti[key]) for key in betti.keys()])
    outBetti.close()
    outSeries = open(betti_dir+'/bettiSeries_{}_{}.txt'.format(ind[0],ind[1]),"w+")
    f="+".join(["{}*t_0^({})*t_1^({})*s_0^({})*s_1^({})".format(betti[key],key[0],key[1],key[2],key[3]) for key in betti.keys()])
    outSeries.write(f)
    outSeries.close()
    print(ind)






    
