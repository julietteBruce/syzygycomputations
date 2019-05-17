#!/usr/bin/python3
import argparse
import subprocess
import tempfile
import os
import subprocess
import os.path

def createMatrices_pq(n,d,b,p,q,out_dir):
   subprocess.check_call(["./build/DirectConstructMatrices","--product",str(n[0]),str(n[1]),str(d[0]),str(d[1]),str(p),out_dir, str(q), str(b[0]), str(b[1])] )


#note binom(d+n,d) is the dimension of S_d, and there is no matrix after that

def createMatrices(n, d, b, pqs, out_dir):
    matDir={};
    for (p,q) in pqs:
        curr_dir=os.path.join(out_dir,"map_{}_{}".format(p,q))
        if not os.path.isdir(curr_dir):
            os.makedirs(curr_dir);
        matDir[(p,q)] = curr_dir
        createMatrices_pq(n,d,b,p,q,curr_dir);
    return matDir

