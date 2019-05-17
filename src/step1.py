import argparse
import os
import os.path
import subprocess


from setup_matrices_prod_Pn import *
from compute_rank_magma import *

argparser = argparse.ArgumentParser();
argparser.add_argument('output_dir')
argparser.add_argument('d1', type = int)
argparser.add_argument('d2', type = int)
argparser.add_argument('b1', type = int)
argparser.add_argument('b2', type = int)

args = argparser.parse_args()

d1=args.d1
d2=args.d2
b1=args.b1
b2=args.b2



out_dir = args.output_dir
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)



M2_input_file = open(os.path.join(out_dir,"computeRR.m2"),'w')
M2_input_file.writelines(['load "relevantRange.m2"\n',
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
M2_input_file.close()


## have relevantRange.m2 input this file. 

m2Path = open(os.path.join(out_dir,"computeRR.m2"));
subprocess.run("M2", stdin=m2Path);
m2Path.close()

## have M2 output "relevantpq.txt"

PnFile=open(os.path.join(out_dir,"relevantpq.txt"),'r')

n=[int(i) for i in (PnFile.readline()).split()];
d=[int(i) for i in (PnFile.readline()).split()];
b=[int(i) for i in (PnFile.readline()).split()];
pqs=[map(int,l.split()) for l in PnFile];




matrix_dir = os.path.join(args.output_dir,"matrices")
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)


matDirs = createMatrices(n, d, b, pqs, matrix_dir)

rankDict={}
for ((p,q),matDir) in matDirs.items():
    rankDict[(p,q)] = call_magma_dir(matDir);

print(rankDict)

