#!/usr/bin/python3
import argparse
import os
import os.path
import cmd
import ast
import traceback
import re
from rank_backends import *

#given a folder with all the matricies in it, produce a folder with the relevant ranks, and then the relevant betti table entries

def binom(n,k):
    val = 1
    for i in range(0,k):
        val*=n-i
    for i in range(1,k+1):
        val/=i
    return int(val)

#the number of derangements of md
def multiplicity(md):
    counts = dict()
    for i in md:
        if i not in counts:
            counts[i] = 0
        counts[i] += 1
    n = len(md)
    mult = 1
    for (i,c) in counts.items():
        mult*=binom(n,c)
        n-=c
    return mult

class BettiTableState(object):
    """
    this represents the state of the betti table computation, in particular,
    the betti table could be partially computed, and this partial state is
    stored on disk

    The structure of the state is a directory with the following layout
    state/
      info.txt - contains one line with N and P and K
      table.txt - a text file representing the betti table.
      matricies/ - contains all the matricies in folders sorted by which map they come fom
        map_p_q/ - the directory containing the matricies for the p,q-th map
      ranks/ - mirrors the format of matricies but with rank data
    
    The program rank_prog should accept two file names as arguments,
    one is for a directory containing all the entries,
    the second is a for a file to output all the ranks into
    """
    def __init__(self,path,rank_backend):
        info_path = os.path.join(path,'info.txt')
        self.__path = path
        self.__matrix_dir = os.path.join(path,'matricies');
        self.__ranks_dir = os.path.join(path,'ranks');
        if not os.path.exists(self.__ranks_dir):
            os.mkdir(self.__ranks_dir)
        if not (os.path.isfile(info_path) and
                os.path.isdir(self.__matrix_dir) and
                os.path.isdir(self.__ranks_dir)):
            raise -1;
        with open(info_path) as info_file:
            info = info_file.readline()
            nstr,dstr,kstr = info.split()
            self.__N = int(nstr)
            self.__D = int(dstr)
            self.__K = int(kstr)
        #Note, table.txt is not loaded since we use that only as output
        self.__rank_backend = rank_backend

        self.__max_q = self.__N
        self.__max_p = binom(self.__D+self.__N,self.__N)

    #deletes all of the relevant ranks
    def clean(self):
        files = os.listdir(self.__ranks_dir);
        for f in files:
            os.remove(os.path.join(self.__ranks_dir,f));
        if os.path.exists(os.path.join(self.__path,'table.txt')):
            os.remove(os.path.join(self.__path,'table.txt'))
    def compute_rank(self,p,q):
        if p<=0 or q<0 or p>self.__max_p or q>self.__max_q:
            return 0
        base_name = "map_{}_{}".format(p,q)
        if not os.path.exists(os.path.join(self.__ranks_dir,base_name + ".txt")):
            self.__rank_backend.compute_ranks(os.path.join(self.__matrix_dir,base_name),
                                              os.path.join(self.__ranks_dir,base_name+".txt"))
        return self.lookup_rank(p,q)
    def lookup_rank(self,p,q):
        if p<=0 or q<0 or p>self.__max_p or q>self.__max_q:
            return 0
        base_name = os.path.join(self.__ranks_dir,"map_{}_{}".format(p,q))
        if os.path.exists(base_name + ".txt"):
            with open(base_name + ".txt") as f:
                rank = 0
                for line in f:
                    match = re.match('(\(.*\)) (\d*)',line)
                    mdstr,rankstr = match.groups()
                    md = ast.literal_eval(mdstr)
                    rank += multiplicity(md)*int(rankstr)
                return rank
        else:
            return None
    def get_betti_entry(self,p,q):
        #the space is is \Wedge^{p}S_d \otimes S_(qd+k)
        #note the size of the basis of S_d is binom(d+n,n)
        #and then the \Wedge is binom(|S_d|,p)
        #and then the \tensor takes the cartesian product of the bases
        basis_size = binom(binom(self.__D+self.__N,self.__N),p)*binom(self.__D*q+self.__K+self.__N,self.__N);
        ker_rank = basis_size - self.compute_rank(p,q)
        img_rank = self.compute_rank(p+1,q-1)
        return ker_rank-img_rank
    #TODO reduce the copy-paste between versions
    def get_betti_table(self):
        return [[self.get_betti_entry(p,q) for p in range(0,self.__max_p+1)] for q in range(0,self.__max_q+1)]
    def lookup_betti_table(self):
        return [[self.lookup_betti_entry(p,q) for p in range(0,self.__max_p+1)] for q in range(0,self.__max_q+1)]
    def lookup_betti_entry(self,p,q):
        basis_size = binom(binom(self.__D+self.__N,self.__N),p)*binom(self.__D*q+self.__N,self.__N);
        outwards_rank = self.lookup_rank(p,q)
        if outwards_rank==None:
            return None
        ker_size = basis_size - outwards_rank
        img_size = self.lookup_rank(p+1,q-1)
        if img_size==None:
            return None
        return ker_size-img_size


class BettiCmd(cmd.Cmd):
    def __init__(self,state):
        super(BettiCmd,self).__init__()
        self.__state = state
        self.prompt = ">"
    def do_print(self,args):
        """print the betti table"""
        argv = args.split()
        if len(argv)==0:
            print(self.__state.lookup_betti_table())
        elif len(argv)==2:
            p,q = map(int,argv)
            print(self.__state.lookup_betti_entry(p,q))
        else:
            print("print expects either zero or two integers as arguments")
        return False
    def do_compute(self,args):
        """compute the betti table"""
        try:
            argv = args.split()
            if len(argv)==0:
                print(self.__state.get_betti_table())
            elif len(argv)==2:
                p,q = map(int,argv)
                print(self.__state.get_betti_entry(p,q))
            else:
                print("compute expects either zero or two integers as arguments")
        except:
            traceback.print_exc()
        return False
    def do_clean(self,args):
        """cleans all cached rank information"""
        self.__state.clean()
        return False
    def do_exit(self,args):
        """quits the program"""
        return True
    def emptyline(self):
        return self.do_help("")

argparser = argparse.ArgumentParser();

argparser.add_argument('path')
argparser.add_argument('--rank-backend')

args = argparser.parse_args()
backend = None
if args.rank_backend=="sage":
    backend = SageBackend()
elif args.rank_backend=="matlab-qr":
    backend = MatlabQRBackend()
elif args.rank_backend=="matlab-lu":
    backend = MatlabLUBackend()
elif args.rank_backend==None:
    print("No backend specified, using sage as a default");
    backend = SageBackend()
else:
    print("unrecognized backend '{}'\n".format(args.rank_backend));
    exit(1)

state = BettiTableState(args.path,backend)

command = BettiCmd(state)
command.cmdloop()
