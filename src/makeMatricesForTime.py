import argparse
import os
import os.path
import subprocess
import re
import time

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
computeRR_file.writelines(['load "relevantRange.m2";\n',
                          "d={"+"{},{}".format(str(d1),str(d2)) + "};\n",
                          "b={"+"{},{}".format(str(b1),str(b2)) + "};\n",
                          "rR=relevantRange(0,b,d);\n", 
                          'g=openOut "{}/relevantpq.txt";\n'.format(out_dir),
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


matrix_dir = os.path.join(args.output_dir,"matrices")



if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)



start = time.time()
matDirs = createMatrices(n, d, b, matricesToCompute, matrix_dir)
stop = time.time()

print(stop-start)

'''
def matrixSize(matFile):
    with open(matFile,"r") as mF:
        sM = mF.readlines()
    entries = list(map(lambda s: list(map(lambda i: int(i), s.split())), sM ))
    rows = set(map(lambda s:s[0], entries))
    columns = set(map(lambda s:s[1], entries))
    return (len(rows),len(columns))


def largestMatrixPQ(p,q):
    pqDir = os.path.join(matrix_dir,"map_{}_{}".format(p,q))
    pqFiles = os.listdir(pqDir)
    if len(pqFiles) == 0:
        return (0,0)
    else:
        rowsCols = list(map(lambda f: matrixSize(os.path.join(pqDir,f)), pqFiles))
        return max(rowsCols, key=lambda rc: rc[0]*rc[1])



numMatrices = sum([ len(os.listdir(os.path.join(matrix_dir,"map_{}_{}".format(p,q)))) for (p,q) in matricesToCompute])

largestMatrix = max([ largestMatrixPQ(p,q) for (p,q) in matricesToCompute], key=lambda rc:rc[0]*rc[1])
print(numMatrices)
print(largestMatrix)


'''
