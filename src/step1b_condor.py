import argparse
import os
import os.path
import subprocess
import glob


from setup_matrices_prod_Pn import *
from compute_rank_magma import *
from ranksToBetti import *


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

rankDict={}
for (p,q) in matricesToCompute:
    rankDict[(p,q)] = {}
    ranks_pq_dir = os.path.join(ranks_dir,"map_{}_{}".format(p,q))
    for mdFile in glob.glob(os.path.join(ranks_pq_dir, "*.ranks")):
        md = filename_to_md(mdFile)
        rankDict[(p,q)][md] =  file_to_rank(mdFile)



betti_dir = os.path.join(args.output_dir,"betti")
if not os.path.isdir(betti_dir):
    os.makedirs(betti_dir)

bettiDict={}
for (p,q) in bettiToCompute:
    bettiPQ=betti_pq(p,q,rankDict)
    bettiDict[(p,q)]=bettiPQ
    with open(os.path.join(betti_dir,'bettiMulti_{}_{}.txt'.format(p,q)),"w+") as outBetti:
        for md in bettiPQ.keys():
            outBetti.write('{} {} {} {} {}\n'.format(md[0],md[1],md[2],md[3],bettiPQ[md]))
    with open(os.path.join(betti_dir,'bettiSeries_{}_{}.txt'.format(p,q)),"w+") as outSeries:
        f="+".join(["{}*t_0^({})*t_1^({})*t_2^({})*t_3^({})".format(bettiPQ[md],md[0],md[1],md[2],md[3]) for md in bettiPQ.keys() if bettiPQ[md]!=0])
        outSeries.write(f)
    with open(os.path.join(betti_dir,'bettiTotal_{}_{}.txt'.format(p,q)),"w+") as outTotal:
        outTotal.write(str(sum([bettiPQ[md] for md in bettiPQ.keys()])))

    


subprocess.run(["tar", "-czf", os.path.join(args.output_dir,"matrices_{}_{}_{}_{}.tar.gz".format(d1,d2,b1,b2)) , matrix_dir])
subprocess.run(["tar", "-czf", os.path.join(args.output_dir,"ranks_{}_{}_{}_{}.tar.gz".format(d1,d2,b1,b2)) , ranks_dir])

subprocess.run(["rm", "-rf", matrix_dir, ranks_dir])
