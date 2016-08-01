#!/usr/bin/python3
import argparse
import os
import os.path

def create_submit(dataDir,p,q):
    template ="""executable = wrapper.sh
output = outdir/single_entry_{1}_{2}.$(CLUSTER).$(PROCESS).out
error = outdir/single_entry_{1}_{2}.$(CLUSTER).$(PROCESS).err
log = single_entry_{1}_{2}.$(CLUSTER).log

### from slide 21 of
##    https://twiki.opensciencegrid.org/twiki/pub/Education/UserSchool16Materials/gthain-2016-HTCondor-3.ppt
request_memory = ifthenelse(MemoryUsage =!= undefined,max({{2048, (MemoryUsage * 3/2)}}), 2048)
periodic_hold = (MemoryUsage >= ((RequestMemory) * 5/4 )) && (JobStatus == 2)
periodic_release = (JobStatus == 5) && ((CurrentTime - EnteredCurrentStatus) > 180) && (NumJobStarts < 5) && (HoldReasonCode =!= 13) && (HoldReasonCode =!= 34)
periodic_remove = (NumJobStarts > 10)
on_exit_remove = (ExitCode =?= 0) && (ExitBySignal =?= false)
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
