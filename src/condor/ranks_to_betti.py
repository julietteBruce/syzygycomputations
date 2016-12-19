#!/usr/bin/python3
import sys
import argparse
import os
import os.path

argparser = argparse.ArgumentParser();
argparser.add_argument('p')
argparser.add_argument('q')

args = argparser.parse_args()

p=int(args.p)
q=int(args.q)

if q!=0:
    print("This code only works for the first row!!!")
    exit(1)

with open("./info.txt") as info:
    (n,d,k) = map(int,info.read().split());

def md_to_filename(md,ext):
    """Takes a multidegree and an extension an creates a file name"""
    return "multidegree_{}{}".format("_".join(map(str,md)),ext)
    
def compute_betti(md,rank):
    matFilename = os.path.join("./matricies/map_{0}_{1}/".format(p,q),md_to_filename(md,".dat"));
    best = 0
    with open (matFilename) as f:
        for line in f:
            (r,c,v) = map(int,line.split())
            if r>best:
                best = r
        md_basis_size = best
        ker_rank = md_basis_size - rank
        return ker_rank;

for line in sys.stdin:
    terms = list(map(int,line.split()))
    md = terms[:-1]
    rank = terms[-1]
    betti = compute_betti(md,rank)
    print("{} {}".format(" ".join(map(str,md)),betti))
