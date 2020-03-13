import argparse
import os
import os.path
import subprocess
import glob
from multiprocessing import Pool
from functools import partial

from setup_matrices_prod_Pn import *
from compute_rank_magma import *
from ranksToBetti import *
#from check_ranks import *

argparser = argparse.ArgumentParser();
argparser.add_argument('output_dir')
argparser.add_argument('d1', type = int)
argparser.add_argument('d2', type = int)
argparser.add_argument('b1', type = int)
argparser.add_argument('b2', type = int)
argparser.add_argument('char', type = int)

args = argparser.parse_args()

d1=args.d1
d2=args.d2
b1=args.b1
b2=args.b2
char = args.char

out_dir = args.output_dir
matrix_dir = os.path.join(args.output_dir,"matrices")
ranks_dir=os.path.join(args.output_dir,"ranks")
condor_dir = os.path.join(args.output_dir,"condor")
submit_dir = os.path.join(condor_dir,"submit")
log_dir = os.path.join(condor_dir, "log")
error_dir = os.path.join(condor_dir, "error")

betti_dir = os.path.join(args.output_dir,"betti")
if not os.path.isdir(betti_dir):
    os.makedirs(betti_dir)


PnFile=open(os.path.join(out_dir,"relevantpq.txt"),'r')

n=[int(i) for i in (PnFile.readline()).split()];
d=[int(i) for i in (PnFile.readline()).split()];
b=[int(i) for i in (PnFile.readline()).split()];
pqs=[tuple(map(int,l.split())) for l in PnFile];

sumPQ = set([p+q for (p,q) in pqs])
pqs_split = [sorted([(p,q) for (p,q) in pqs if p+q==r], key = lambda tup: tup[1]) for r in sumPQ]

bettiToCompute = []
for pql in pqs_split:
    if len(pql) == 1:
        continue
    else:
        bettiToCompute += pql[:-1] 

matricesToCompute = bettiToCompute.copy()
for (p,q) in bettiToCompute:
    if (p+1,q-1) not in bettiToCompute:
        matricesToCompute.append((p+1,q-1))


def filename_to_md(name):
    """Take a filename of the format *multidegree_* and reads off the multidegree"""
    match = re.search('multidegree((?:_\d+)*)',name)
    return tuple(map(int,re.findall("\d+",match.group(1))))
        
def file_to_rank(name):
    with open(name,'r') as f:
        fLines = f.readlines()
        (rank, cols) = tuple(map(int,fLines))
    return (rank, cols)

def multidegrees(p,q):
    matrix_pq_dir = os.path.join(matrix_dir, "map_{}_{}".format(p,q))
    return [filename_to_md(mdFile)  for mdFile in glob.glob(os.path.join(matrix_pq_dir, "*.dat"))]


def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def remove_empty_rank_files(ranks_pq_dir):
    for mdFile in glob.glob(os.path.join(ranks_pq_dir, "*.ranks")):
        if os.path.getsize(mdFile)==0:
            os.remove(mdFile)

    

def missing_ranks(p,q):
    ranks_pq_dir = os.path.join(ranks_dir, "map_{}_{}".format(p,q))
    matrices_md = set(multidegrees(p,q))
    ranks_md = set()
    for mdFile in glob.glob(os.path.join(ranks_pq_dir, "*.ranks")):
        if os.path.getsize(mdFile)==0:
            os.remove(mdFile)
            continue
        with open(mdFile,'r') as f:
            flines = f.readlines()
        if len(flines) != 2:
            os.remove(mdFile)
            continue
        else:
            if not RepresentsInt(flines[0]) or not RepresentsInt(flines[1]):
                os.remove(mdFile)
                continue
            else:
                ranks_md.add(filename_to_md(mdFile))
    return matrices_md - ranks_md



def make_rank_dict(matToComp):
    rankDict={}
    
    for (p,q) in matToComp:
        
        rankDict[(p,q)] = {}
        ranks_pq_dir = os.path.join(ranks_dir,"map_{}_{}".format(p,q))
        
        for mdFile in glob.glob(os.path.join(ranks_pq_dir, "*.ranks")):
            md = filename_to_md(mdFile)
            rankDict[(p,q)][md] =  file_to_rank(mdFile)

        print(multidegrees(p,q))
        missing = set(multidegrees(p,q)) - set(rankDict[(p,q)].keys())
        print(missing)
        for md in missing:
            rankDict[(p,q)][md] = "infinity";
    return rankDict

    


def make_betti(matToComp, bettiToComp):
    rankDict = make_rank_dict(matToComp)
    bettiDict={}
    for (p,q) in bettiToComp:
        bettiPQ=betti_pq(p,q,rankDict)
        bettiDict[(p,q)]=bettiPQ
        with open(os.path.join(betti_dir,'bettiMulti_{}_{}.txt'.format(p,q)),"w+") as outBetti:
            for md in bettiPQ.keys():
                outBetti.write('{} {} {} {} {}\n'.format(md[0],md[1],md[2],md[3],bettiPQ[md]))
        with open(os.path.join(betti_dir,'bettiSeries_{}_{}.txt'.format(p,q)),"w+") as outSeries:
            f="+".join(["{}*t_0^({})*t_1^({})*t_2^({})*t_3^({})".format(bettiPQ[md],md[0],md[1],md[2],md[3]) for md in bettiPQ.keys() if bettiPQ[md]!=0])
            outSeries.write(f)
        with open(os.path.join(betti_dir,'bettiTotal_{}_{}.txt'.format(p,q)),"w+") as outTotal:
            if all([RepresentsInt(bettiPQ[md]) for md in bettiPQ.keys()]):
                outTotal.write(str(sum([bettiPQ[md] for md in bettiPQ.keys()])))
            else:
                outTotal.write("infinity")
    return bettiDict
    

mb = make_betti(matricesToCompute, bettiToCompute)

