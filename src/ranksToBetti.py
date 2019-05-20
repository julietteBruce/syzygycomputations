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

HARDCODE1=ranks_dir
HARDCODE2=betti_dir


#HARDCODE1 = 'mgRankData'
#HARDCODE2 = 'mgBettiData'

onlyfiles = glob.glob(os.path.join(HARDCODE1,"*.txt"))
onlyfiles_base = [os.path.basename(f) for f in onlyfiles]
indices=[f.split('_')[1:-1] for f in onlyfiles_base]
indices=[list(map(int,s)) for s in indices]


#onlyfiles = [f for f in listdir(HARDCODE1) if isfile(join(HARDCODE1, f))]
#del onlyfiles[0]
#indices = [re.split("_", s[4:-4]) for s in onlyfiles]
#indices = [[int(s[0]),int(s[1])] for s in indices]

for ind in indices:
    d = {}
    with open(HARDCODE1+'/ranks_{}_{}.txt'.format(ind[0],ind[1])) as f:
        for line in f:
            (key, val) = line.split(')')
            d[eval(key+')')] = map(int, val[1:-1].split())
    f.close
    if ind[1] != 0:
        E = {}
        with open(HARDCODE1+'/ranks_{}_{}.txt'.format(ind[0]+1,ind[1]-1)) as f:
            for line in f:
                (key, val) = line.split(')')
                E[eval(key+')')] = map(int, val[1:-1].split())
        f.close
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
    outBetti = open(HARDCODE2+'/betti_{}_{}.txt'.format(ind[0],ind[1]),"w+")
    for key in betti.keys():
        if key == betti.keys()[-1]:
            f = str(key)+' '+str(betti[key])
            outBetti.write(f)
        else: 
            f = str(key)+' '+str(betti[key])+'\n'
            outBetti.write(f)
    outBetti.close()
    outSeries = open(HARDCODE2+'/bettiSeries_{}_{}.txt'.format(ind[0],ind[1]),"w+")
    for key in betti.keys():
        if key == betti.keys()[-1]:
            fx = 't_0^({})+t_1^({})'.format(key[0],key[1])
            fy = 's_0^({})*s_1^({})'.format(key[2],key[3])
            f  = str(betti[key])+'*'+fx+'*'+fy
            outSeries.write(f)
        else:
            fx = 't_0^({})*t_1^({})'.format(key[0],key[1])
            fy = 's_0^({})*s_1^({})'.format(key[2],key[3])
            f  = str(betti[key])+'*'+fx+'*'+fy+'+'
            outSeries.write(f)
    outSeries.close()
    print(ind)
