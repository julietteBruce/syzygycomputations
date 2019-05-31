import argparse
import os
import os.path
import subprocess



from setup_matrices_prod_Pn import *
from compute_rank_magma import *
from ranksToBetti import *

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

PnFile=open(os.path.join(out_dir,"relevantpq.txt"),'r')

n=[int(i) for i in (PnFile.readline()).split()];
d=[int(i) for i in (PnFile.readline()).split()];
b=[int(i) for i in (PnFile.readline()).split()];
pqs=[tuple(map(int,l.split())) for l in PnFile];
for pq in pqs:
    if pq[1]==1 and (pq[0]+1,pq[1]-1) not in pqs:
        pqs.append((pq[0]+1, pq[1]-1))

matrix_dir = os.path.join(args.output_dir,"matrices")
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)

ranks_dir=os.path.join(args.output_dir,"ranks")
if not os.path.isdir(ranks_dir):
    os.makedirs(ranks_dir)

matDirs = createMatrices(n, d, b, pqs, matrix_dir)

rankDict={}
for ((p,q),matDir) in matDirs.items():
    if q==2: ## since we only need to know (p,1) or (p-1,2) for the b=0 case. Fix for b!=0 case
        continue
    ranks_p_q_file=open(os.path.join(ranks_dir,"ranks_{}_{}.txt".format(p,q)),'w')
    rankDict[(p,q)] = call_magma_dir(matDir);
    ranks_p_q_file.writelines([' '.join(list(map(str,md)) + [str(rankDict[(p,q)][md][0]), str(rankDict[(p,q)][md][1])+'\n']) for md in rankDict[(p,q)].keys()])



betti_dir = os.path.join(args.output_dir,"betti")
if not os.path.isdir(betti_dir):
    os.makedirs(betti_dir)


bettiDict={}
for pq in pqs:
    if pq[1]==2: ## since we only need to know (p,1) or (p-1,2) for the b=0 case. Fix for b!=0 case
        continue
    outBetti = open(os.path.join(betti_dir,'betti_{}_{}.txt'.format(pq[0],pq[1])),"w+")
    bettiPQ=betti_pq(pq[0],pq[1],rankDict)
    bettiDict[pq]=bettiPQ
    for md in bettiPQ.keys():
        outBetti.write('{} {} {} {} {}\n'.format(md[0],md[1],md[2],md[3],bettiPQ[md]))
        if md[0]<md[1]:
            outBetti.write('{} {} {} {} {}\n'.format(md[1],md[0],md[2],md[3],bettiPQ[md]))
            if md[2]<md[3]:
                outBetti.write('{} {} {} {} {}\n'.format(md[1],md[0],md[3],md[2],bettiPQ[md]))
        if md[2]<md[3]:
            outBetti.write('{} {} {} {} {}\n'.format(md[0],md[1],md[3],md[2],bettiPQ[md]))
    outBetti.close()
    outSeries=open(os.path.join(betti_dir,'bettiSeries_{}_{}.txt'.format(pq[0],pq[1])),"w+") # no longer gives correct answer, to fix
    f="+".join(["{}*t_0^({})*t_1^({})*s_0^({})*s_1^({})".format(bettiPQ[md],md[0],md[1],md[2],md[3]) for md in bettiPQ.keys()])
    outSeries.write(f)
    outSeries.close()

    


#betti_ret = subprocess.run(["python3", "src/ranksToBetti.py" , ranks_dir, betti_dir])


