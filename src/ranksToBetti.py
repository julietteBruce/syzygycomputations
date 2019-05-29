import numpy as np
import sys
import re
import itertools
from os import listdir
from os.path import isfile, join
import os
import argparse
import glob


'''
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

def rankDict_pq(p,q,):
    rankDictpq={}
    f=open(os.path.join(ranks_dir,'ranks_{}_{}.txt'.format(p,q),'r'))
    f_lines=f.readlines()
    for line in f_lines:
        dat = list(map(int, line.split(' ')))
        rankDictpq[tuple(dat[:4])] = tuple(dat[4:])
    f.close()
    return rankDictpq

def rankDict(pq_list):
    return {pq:rankDict_pq(pq[0],pq[1]) for pq in pq_list}

'''


def betti_pq(p,q,rank_dict):
    bettiPQ={}
    for md in rank_dict[(p,q)].keys():
        if q-1<0:
            bettiPQ[md] = rank_dict[(p,q)][md][1] - rank_dict[(p,q)][md][0]
        else:
            if md not in rank_dict[(p+1,q-1)].keys():
                bettiPQ[md] =  rank_dict[(p,q)][md][1] - rank_dict[(p,q)][md][0]
            else:
                bettiPQ[md] = rank_dict[(p,q)][md][1] - rank_dict[(p,q)][md][0] - rank_dict[(p+1,q-1)][md][0]
    return bettiPQ



'''
def betti(pq_list):
    bettiDict={}
    for pq in pq_list:
        outBetti = open(os.path.join(betti_dir,'betti_{}_{}.txt'.format(pq[0],pq[1])),"w+")
        bettiPQ=betti_pq(pq[0],pq[1],pq_list)
        bettiDict[pq]=bettiPQ
        for md in bettiPQ.keys():
           outBetti.write('{} {}\n'.format(md,bettiPQ[md]))
        outBetti.close()
        outSeries=open(os.path.join(betti_dir,'bettiSeries_{}_{}.txt'.format(pq[0],pq[1])),"w+")
        f="+".join(["{}*t_0^({})*t_1^({})*s_0^({})*s_1^({})".format(betti[pq][md],md[0],md[1],md[2],md[3]) for md in betti[pq].keys()])
        outSeries.write(f)
        outSeries.close()
        return bettiDict

           
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
'''
