import sys
import ast
import glob
import os.path

def load_matrix(file):
    entries = dict();
    for line in file:
        x,y,val = map(int,line.split())
        entries[(x,y)] = val
    m = matrix(ZZ,entries,sparse=True);
    rows = m.columns();
    def nonzero(v):
        return v.parent().zero()!=v;
    return matrix(ZZ,filter(nonzero,rows),sparse=False)
    


for f in glob.glob(os.path.join(sys.argv[1],"*.dat")):
    base_name = os.path.splitext(os.path.basename(f))[0]
    parts = base_name.split('_')
    md = tuple(map(int,parts[1:]))
    with open(f) as file:
        print("{} {}".format(md,load_matrix(file).rank()))
