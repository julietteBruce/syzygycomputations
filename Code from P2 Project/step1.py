#!/usr/bin/python3
import sys
import argparse
import os
import os.path
import re
import itertools

argparser = argparse.ArgumentParser();
argparser.add_argument('n',type=int);
argparser.add_argument('e',type=int);
argparser.add_argument('b',type=int);
argparser.add_argument('p',type=int);
argparser.add_argument('q',type=int);

args = argparser.parse_args();

p = args.p;
q = args.q;


def filename(base,prefix,p,q):
    return "{}{}_{}/{}{}_{}.txt".format(base,args.e,args.b,prefix,p,q);

outfile = filename("./mgBettiData/P{}/".format(args.n),"",p,q);

inDirectory = "./mgRankData/P{}/".format(args.n);

ker_rank_file = filename(inDirectory,"ranks_",args.p,args.q) 
if q!=0:
    img_rank_file = ker_rank_file = filename(inDirectory,"ranks_",p+1,q-1);
dimension_file = filename(inDirectory,"dimensions_",args.p,args.q);

dimensions = dict({});

with open(dimension_file) as f:
    for line in f:
        match = re.search('(\(.*\)) (\d+)',line)
        md = tuple(map(int,re.findall("\d+",match.group(1))))
        val = int(match.group(2))
        dimensions[md] = val;

def read_rank_file(file):
    ranks = dict({});
    with open(ker_rank_file) as f:
        for line in f:
            terms = tuple(map(int,line.split()));
            md = terms[:-1];
            rank = terms[-1];
            ranks[md] = rank;
    return ranks;

ker_ranks = read_rank_file(ker_rank_file);
if q!=0:
    img_ranks = read_rank_file(img_rank_file);


try:
    os.mkdir(os.path.dirname(outfile));
except:
    None
    
def all_permutations(lst):
    return set(itertools.permutations(lst));
    
with open(outfile,'w') as file:
    for md in ker_ranks:
        betti = dimensions[md]-ker_ranks[md] - (img_ranks[md] if q!=0 else 0);
        for real_md in all_permutations(md):
            file.write("{},{}\n".format(",".join(map(str,real_md)),betti));
    
