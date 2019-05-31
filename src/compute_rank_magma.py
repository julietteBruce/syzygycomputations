# python file to make magma files to compute ranks

import re
import glob
import os.path
import subprocess



def filename_to_md(name):
    """Take a filename of the format *multidegree_* and reads off the multidegree"""
    match = re.search('multidegree((?:_\d+)*)',name)
    return tuple(map(int,re.findall("\d+",match.group(1))))



def call_magma_file(matrixFile):
    magmaPath = os.path.join(os.path.dirname(__file__),"ranks.magma");
    ret = subprocess.run(["magma","-b","file:=" + matrixFile,magmaPath],stdout = subprocess.PIPE,check=True)
    if ret.returncode==0 and len(ret.stdout)!=0:
        #print(ret.stdout)
        return list(map(int, ret.stdout.decode().split('\n')[:-1]))
    else:
        raise Exception("Couldn't run magma");

def call_magma_dir(matrixDir):
    rankDict={}
    for fileName in glob.glob(os.path.join(matrixDir,"*.dat")):
        md=filename_to_md(fileName)
        if md[0]>md[1] or md[2]>md[3]: ## remove redundant multidegrees - will be unnecessary when we remove redundant multidegrees from matrices construction
            continue
        rankDict[filename_to_md(fileName)] = call_magma_file(fileName)
    return rankDict
    

    
