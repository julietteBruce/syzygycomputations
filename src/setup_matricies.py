import argparse
import subprocess
import tempfile
import os
import os.path

def binom(n,k):
    val = 1
    for i in range(0,k):
        val*=n-i
    for i in range(1,k+1):
        val/=i
    return int(val)

argparser = argparse.ArgumentParser();
argparser.add_argument('output_dir')
argparser.add_argument('n')
argparser.add_argument('d')

args = argparser.parse_args()

d=int(args.d)
n=int(args.n)

matrix_dir = os.path.join(args.output_dir,"matricies")
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)

#note binom(d+n,d) is the dimension of S_d, and there is no matrix after that
for p in range(1, binom(d+n,n)+1):
    with tempfile.NamedTemporaryFile() as temp:
        subprocess.check_call(["./build/ConstructMatricies",args.n,args.d,str(p)],stdout=temp)
        for q in range(0,n+1):
            curr_dir=os.path.join(matrix_dir,"map_{}_{}".format(p,q))
            if not os.path.isdir(curr_dir):
                os.makedirs(curr_dir);
            subprocess.check_call(["./build/SliceMatrix",temp.name,str(q),curr_dir])

with open(os.path.join(args.output_dir,"info.txt"),"w") as info_file:
    info_file.write("{} {}".format(n,d));

