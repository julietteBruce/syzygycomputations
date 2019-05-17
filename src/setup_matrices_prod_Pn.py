#!/usr/bin/python3
import argparse
import subprocess
import tempfile
import os
import os.path
import itertools
import time


argparser = argparse.ArgumentParser();
argparser.add_argument('output_dir')
argparser.add_argument('Pns')

args = argparser.parse_args()

PnFile=open(args.Pns,'r')

ns=[int(i) for i in (PnFile.readline()).split()] 
ds=[int(i) for i in (PnFile.readline()).split()]
bs=[int(i) for i in (PnFile.readline()).split()]

pqs=[map(int,l.split()) for l in PnFile]




#entries_set = set(map(tuple,args.entries))
#if len(args.entries)!=0:
#    p_set = {p for (p,q) in args.entries if p<=binom(d+n,n)}
#    p_set |= {p+1 for p in p_set if p+1<=binom(d+n,n)}
#else:
#    p_set = set(range(1, binom(d+n,n)+1));




matrix_dir = os.path.join(args.output_dir,"matrices")
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)

log_dir = os.path.join(args.output_dir,"logs/outdir");
if not os.path.isdir(log_dir):
    os.makedirs(log_dir)



def createMatrices(p,q,bs,out_dir):
   subprocess.check_call(["./build/DirectConstructMatrices","--product",str(ns[0]),str(ns[1]),str(ds[0]),str(ds[1]),str(p),out_dir, str(q), str(bs[0]), str(bs[1])] )


#note binom(d+n,d) is the dimension of S_d, and there is no matrix after that
for (p,q) in pqs:
    curr_dir=os.path.join(matrix_dir,"map_{}_{}".format(p,q))
    if not os.path.isdir(curr_dir):
        os.makedirs(curr_dir);
    createMatrices(p,q,bs,curr_dir);

