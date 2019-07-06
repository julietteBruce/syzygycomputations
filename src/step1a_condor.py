import argparse
import os
import os.path
import subprocess
import re


from setup_matrices_prod_Pn import *
from compute_rank_magma import *
from ranksToBetti import *

argparser = argparse.ArgumentParser();
argparser.add_argument('output_dir')
argparser.add_argument('--tmp', help="store matrix directory in /tmp/dcorey")
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
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

computeRR_file = open(os.path.join(out_dir,"computeRR.m2"),'w')
computeRR_file.writelines(['load "src/relevantRange.m2"\n',
                          "d={"+"{},{}".format(str(d1),str(d2)) + "}\n",
                          "b={"+"{},{}".format(str(b1),str(b2)) + "}\n",
                          "rR=relevantRange(0,b,d)\n", 
                          'g=openOut "{}/relevantpq.txt"\n'.format(out_dir),
                          'g<< "1 1" << endl;\n',
                          'g<< concatenate(toString(d#0), " ", toString(d#1)) << endl;\n',
                          'g<< concatenate(toString(b#0), " ", toString(b#1)) << endl;\n',
                          'for i from 0 to #rR-1 do (\n',
                          '    g<<concatenate(toString(rR#i#0), " ", toString(rR#i#1)) << endl;\n',
                          '    );\n',
                          'g<<close;'])    
computeRR_file.close()


## have relevantRange.m2 input this file. 

m2Path = open(os.path.join(out_dir,"computeRR.m2"));
subprocess.run("M2", stdin=m2Path);
m2Path.close()

## have M2 output "relevantpq.txt"

with open(os.path.join(out_dir,"relevantpq.txt"),'r') as PnFile:
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

if args.tmp:
    matrix_dir = '/tmp/dcorey/matrices_{}_{}_{}_{}'.format(d1,d2,b1,b2)
else:
    matrix_dir = os.path.join(args.output_dir,"matrices")

if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)

ranks_dir=os.path.join(args.output_dir,"ranks")
if not os.path.isdir(ranks_dir):
    os.makedirs(ranks_dir)

matDirs = createMatrices(n, d, b, matricesToCompute, matrix_dir)


# ########## Condor option ####################


condor_dir = os.path.join(args.output_dir,"condor")
if not os.path.isdir(condor_dir):
    os.makedirs(condor_dir)
submit_dir = os.path.join(condor_dir,"submit")
if not os.path.isdir(submit_dir):
    os.makedirs(submit_dir)
log_dir = os.path.join(condor_dir, "log")
if not os.path.isdir(log_dir):
    os.makedirs(log_dir)
error_dir = os.path.join(condor_dir, "error")
if not os.path.isdir(error_dir):
    os.makedirs(error_dir)

    
def ranks_condor(p,q):
    ranks_pq_dir = os.path.join(ranks_dir,"map_{}_{}".format(p,q))
    if not os.path.isdir(ranks_pq_dir):
        os.makedirs(ranks_pq_dir)
    with open(os.path.join(submit_dir, "map_{}_{}.submit".format(p,q)),'w') as submitFile:
        submitFile.write("\n".join([
            "universe = vanilla",
            "REQUIREMENTS = IsMagma==True",
            "magmaRanksScript = src/ranks.magma",
            "p = {}".format(char),
            "executable = condor_Math/bin/runRanks.magma.sh",
            "args = $(infile) $(magmaRanksScript) $(p)",
            "outfileTmp = $Fn(infile)",
            "outfile = $Fn(outfileTmp)",
            "output = {}/$(outfile).$(CLUSTER).$(PROCESS).ranks".format(ranks_pq_dir),
            "error = {}/$(outfile).$(CLUSTER).$(PROCESS).err".format(error_dir),
            "log = {}/map_{}_{}.$(CLUSTER).log".format(log_dir,p,q),
            "request_memory = 2G",
            "request_cpus = 1",
            "queue infile matching files {}/map_{}_{}/*.dat".format(matrix_dir,p,q)
        ]))

    subprocess.run(["condor_submit", os.path.join(submit_dir, "map_{}_{}.submit".format(p,q))])

for (p,q) in matricesToCompute:
    ranks_condor(p,q)

# subprocess.run(["rm", os.path.join(out_dir,"computeRR.m2")])
