import numpy as np
import sys
import re
import itertools
from os import listdir
from os.path import isfile, join

HARDCODE1 = 'mgRankData'
HARDCODE2 = 'mgBettiData'

onlyfiles = [f for f in listdir(HARDCODE1) if isfile(join(HARDCODE1, f))]
del onlyfiles[0]
indices = [re.split("_", s[4:-4]) for s in onlyfiles]
indices = [[int(s[0]),int(s[1])] for s in indices]

for ind in indices:
    d = {}
    with open(HARDCODE1+'/out_'+str(ind[0])+'_'+str(ind[1])+'.txt') as f:
        for line in f:
            (key, val) = line.split(')')
            d[eval(key+')')] = map(int, val[1:-1].split())
    f.close
    if ind[1] != 0:
        E = {}
        with open(HARDCODE1+'/out_'+str(ind[0]+1)+'_'+str(ind[1]-1)+'.txt') as f:
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
    outBetti = open(HARDCODE2+'/betti_'+str(ind[0])+'_'+str(ind[1])+'.txt',"w+")
    for key in betti.keys():
        if key == betti.keys()[-1]:
            f = str(key)+' '+str(betti[key])
            outBetti.write(f)
        else: 
            f = str(key)+' '+str(betti[key])+'\n'
            outBetti.write(f)
    outBetti.close()
    outSeries = open(HARDCODE2+'/bettiSeries_'+str(ind[0])+'_'+str(ind[1])+'.txt',"w+")
    for key in betti.keys():
        if key == betti.keys()[-1]:
            fx = 't_0^(' + str(key[0])+')'+'*t_1^(' + str(key[1])+')'
            fy = 's_0^(' + str(key[2])+')'+'*s_1^(' + str(key[3])+')'
            f  = str(betti[key])+'*'+fx+'*'+fy
            outSeries.write(f)
        else:
            fx = 't_0^(' + str(key[0])+')'+'*t_1^(' + str(key[1])+')'
            fy = 's_0^(' + str(key[2])+')'+'*s_1^(' + str(key[3])+')'
            f  = str(betti[key])+'*'+fx+'*'+fy+'+'
            outSeries.write(f)
    outSeries.close()
    print(ind)
