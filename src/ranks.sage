import sys
import ast

def load_matrix(file):
    entries = dict();
    for line in file:
        x,y,val = map(int,line.split())
        entries[(x,y)] = val
    return matrix(ZZ,entries,sparse=False);

for f in os.listdir(sys.argv[1]):
    base_name = f.split('.')[0]
    parts = base_name.split('_')
    md = tuple(map(int,parts[1:]))
    with open(os.path.join(sys.argv[1],f)) as file:
        print("{} {}".format(md,load_matrix(file).rank()))
