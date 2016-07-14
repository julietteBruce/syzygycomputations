#!/usr/bin/python3
import argparse
import os
import os.path

def create_submit(dataDir,p,q):
    template ="""executable = wrapper.sh
output = outdir/single_entry_{1}_{2}.$(CLUSTER).$(PROCESS).out
error = outdir/single_entry_{1}_{2}.$(CLUSTER).$(PROCESS).err
log = single_entry_{1}_{2}.$(CLUSTER).log

universe = vanilla

arguments=$(infile) {0}/out_{1}_{2}/

queue infile matching files {0}/matricies/map_{1}_{2}/*.dat"""
    return template.format(dataDir,p,q);

argparser = argparse.ArgumentParser();
argparser.add_argument('dataDir')
argparser.add_argument('p',type=int)
argparser.add_argument('q',type=int)
args = argparser.parse_args()

p=args.p
q=args.q
out_dir = "{0}/out_{1}_{2}/".format(args.dataDir,p,q)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir);


with open("{0}/single_entry_{1}_{2}.submit".format(args.dataDir,p,q),'w') as f:
    f.write(create_submit(args.dataDir,p,q));
