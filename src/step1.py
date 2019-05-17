
import argparse
import os
import os.path


from setup_matrices_prod_Pn import *
from compute_rank_magma import *
# import compute ranks functions
# magma -b 'file:=/fac/dcorey/Documents/syzygycomputations/test-2-out/matrices/map_1_1/multidegree_0_4_0_6.dat' ranks.magma


# parse input data
# feed to setup matrices -give a list of file names where all matrices (as dictionary)
# pass list of file names to compute ranks - want dictionary of multidegrees to ranks
# output multidegrees somewhere? (one file per p,q? probably not one per multidegree)



argparser = argparse.ArgumentParser();
argparser.add_argument('output_dir')
argparser.add_argument('Pns')

args = argparser.parse_args()

PnFile=open(args.Pns,'r')

n=[int(i) for i in (PnFile.readline()).split()] 
d=[int(i) for i in (PnFile.readline()).split()]
b=[int(i) for i in (PnFile.readline()).split()]
pqs=[map(int,l.split()) for l in PnFile]




matrix_dir = os.path.join(args.output_dir,"matrices")
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)


matDirs = createMatrices(n, d, b, pqs, matrix_dir)

rankDict={}
for ((p,q),matDir) in matDirs.items():
    rankDict[(p,q)] = call_magma_dir(matDir);

print(rankDict)
