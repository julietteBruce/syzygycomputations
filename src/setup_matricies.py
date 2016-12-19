#!/usr/bin/python3
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
argparser.add_argument('n',type = int)
argparser.add_argument('d',type = int)
argparser.add_argument('-k',type = int,default=0)
argparser.add_argument('--entry',dest='entries',help='which entries to construct matricies for',nargs=2,action='append',type=int,default=[])
argparser.add_argument('--EEL-test',dest='eel',help='construct only those matricies needed for EEL weight tests',action="store_true")

args = argparser.parse_args()

d=args.d
n=args.n
k=args.k

entries_set = set(map(tuple,args.entries))
if len(args.entries)!=0:
    p_set = {p for (p,q) in args.entries if p<=binom(d+n,n)}
    p_set |= {p+1 for p in p_set if p+1<=binom(d+n,n)}
else:
    p_set = set(range(1, binom(d+n,n)+1));

matrix_dir = os.path.join(args.output_dir,"matricies")
if not os.path.isdir(matrix_dir):
    os.makedirs(matrix_dir)

print args.eel
if not args.eel:
    def do_slicing(mat,p,q,k,out_dir):
        subprocess.check_call(["./build/SliceMatrix",mat,str(q),str(k),out_dir])
else:
    def do_slicing(mat,p,q,k,out_dir):
        mdListFilename = os.path.join(out_dir,"mdList.txt");
        subprocess.Popen(["M2","--script","./PrintEELConj.m2",str(d),str(p),str(q),str(k),os.path.join("../../",mdListFilename)],cwd="./src/M2Scripts/").wait();
        subprocess.check_call(["./build/SliceMatrix",mat,str(q),str(k),out_dir,mdListFilename]);

#note binom(d+n,d) is the dimension of S_d, and there is no matrix after that
for p in p_set :
    with tempfile.NamedTemporaryFile() as temp:
        subprocess.check_call(["./build/ConstructMatricies",str(args.n),str(args.d),str(p)],stdout=temp)
        for q in range(0,n+1):
            if ((p,q) in entries_set) or ((p-1,q+1) in entries_set) or (len(entries_set)==0):
                curr_dir=os.path.join(matrix_dir,"map_{}_{}".format(p,q))
                if not os.path.isdir(curr_dir):
                    os.makedirs(curr_dir);
                    do_slicing(temp.name,p,q,k,curr_dir);

with open(os.path.join(args.output_dir,"info.txt"),"w") as info_file:
    info_file.write("{} {} {}".format(n,d,k));

