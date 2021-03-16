import sys
import re
import itertools
from os import listdir
from os.path import isfile, join
import os
import argparse
import glob


def filename_to_indices(name):
    """Take a filename of the format *ranks/ranks_* and reads off the multidegree"""
    match = re.search('bettiF0((?:_\d+)*)',name)
    return tuple(map(int,re.findall("\d+",match.group(1))))

def dFirst(bd):
    return tuple([bd[2],bd[3],bd[0],bd[1]])

onlyfiles = glob.glob(os.path.join("HirzebruchSyzygies","*"))
bds = list(set([filename_to_indices(f) for f in onlyfiles]))
bds.sort(key=dFirst)

with open("dataRange.m2","w") as dataRangeFile:
    dataRangeFile.write("{\n")
    dataRangeFile.write(",\n".join(["{{"+"{},{}".format(bd[0],bd[1])+"},{"+"{},{}".format(bd[2],bd[3])+"}}" for bd in bds]))
    dataRangeFile.write("\n}")
