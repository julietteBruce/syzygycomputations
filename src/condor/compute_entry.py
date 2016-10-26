#!/usr/bin/python3
import argparse
import os
import os.path
import glob
import re
import sys
import string

def binom(n,k):
    val = 1
    for i in range(0,k):
        val*=n-i
    for i in range(1,k+1):
        val/=i
    return int(val)

#the number of derangements of md
def multiplicity(md):
    counts = dict()
    for i in md:
        if i not in counts:
            counts[i] = 0
        counts[i] += 1
    n = len(md)
    mult = 1
    for (i,c) in counts.items():
        mult*=binom(n,c)
        n-=c
    return mult

def filename_to_md(name):
    """Take a filename of the format *multidegree_* and reads off the multidegree"""
    match = re.search('multidegree((?:_\d+)*)',name)
    return tuple(map(int,re.findall("\d+",match.group(1))))

def md_to_filename(md,ext):
    """Takes a multidegree and an extension an creates a file name"""
    return "multidegree_{}{}".format("_".join(map(str,md)),ext)

def get_md_rank(p,q,md):
    fname = os.path.join("./out_{0}_{1}/".format(p,q),md_to_filename(sorted(md),".txt"))
    try:
        with open(fname,'r') as f:
            rnk = int(f.read())
            return rnk;
    except IOError:
        outfileName = os.path.join("./matricies/map_{0}_{1}/".format(p,q),md_to_filename(sorted(md),".dat"));
        #if the matrix exists, return Nan otherwise 0
        if os.path.isfile(outfileName):
            return float('NaN');
        else:
            return 0;
    except ValueError:
        return float('NaN');

def get_rank(p,q):
    tot = 0
    isValid = True
    for fname in glob.glob("./out_{0}_{1}/multidegree_*.txt".format(p,q)):
        md = filename_to_md(fname)
        with open(fname,'r') as f:
            try:
                rnk = int(f.read())
                #print(md)
                #print(rnk)
                tot += rnk*multiplicity(md)

            except ValueError:
                isValid = False
                #print(md)
                #print(None)
    if isValid:
        #print(tot)
        return tot
    else:
        #print(None)
        return float('nan');


argparser = argparse.ArgumentParser();

argparser.add_argument('p')
argparser.add_argument('q')
argparser.add_argument('--md')
args = argparser.parse_args()

p=int(args.p)
q=int(args.q)

#read in n d and k from the info file
with open("./info.txt") as info:
    (n,d,k) = map(int,info.read().split());

if(args.md):
    string_md = tuple(args.md.split(','))
    #TODO temporary hack to get some numbers, we need to store the size of matricies
    fname_glob = os.path.join("./matricies/map_{0}_{1}/".format(p,q),md_to_filename(string_md,".dat"))
    for fname in glob.glob(fname_glob):
        md = filename_to_md(fname)
        best = 0
        with open(fname) as f:
            for line in f:
                (r,c,v) = map(int,line.split())
                if r>best:
                    best = r
        md_basis_size = best
        ker_rank = md_basis_size - get_md_rank(p,q,md)
        img_rank = get_md_rank(p+1,q-1,md)
        print("{} {}".format(" ".join(map(str,md)),ker_rank-img_rank))
        sys.stdout.flush()
else:
    basis_size = binom(binom(d+n,n),p)*binom(d*q+k+n,n);
    ker_rank = basis_size - get_rank(p,q)
    img_rank = get_rank(p+1,q-1)
    print(ker_rank-img_rank)
