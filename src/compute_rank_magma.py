#!/usr/bin/python3


# python file to make magma files to compute ranks
import argparse
import sys
import ast
import glob
import os.path


# argparser = argparse.ArgumentParser();
# argparser.add_argument('map_dir')
# argparser.add_argument('rank_dir')
# argparser.add_argument('-field', type=str, default="RationalField()")

# args = argparser.parse_args()


def file_to_magma_sparse(matrixFile,field):
    delta=open(matrixFile,'r')
    entries = [list(map(int,entry.split())) for entry in delta.readlines()]
    rows=max([entry[0] for entry in entries])
    cols=max([entry[1] for entry in entries])
    magma_entries = ','.join(["<{},{},{}>".format(e[0],e[1],e[2]) for e in entries])
    delta.close()
    return "SparseMatrix({},{},{},[{}])".format(field,rows,cols,magma_entries)


def make_magma_file_one_matrix(matrixFile, outFile, field):
    rank_delta = open(outFile, 'w')
    rank_delta.writelines(["d:={};\n".format(file_to_magma_sparse(matrixFile,field)),"Rank(d);"])
    delta.close()
    rank_delta.close()


def make_magma_file_one_map(matrixFiles, outFile, field):
    rank_deltas = open(outFile, 'w')
    sparse_matrices = [file_to_magma_sparse(matrixFile,field) for matrixFile in matrixFiles]
    rank_deltas.write("ds:=[\n" + ",\n".join(sparse_matrices) + "];\n")
    rank_deltas.write("for M in ds do\n"+"print Rank(M);\n"+"end for;")
    rank_deltas.close()


def make_magma_file_from_dir(matrixDir,outFile,field):
    files=[f for f in glob.glob(matrixDir + "/*.dat")]
    make_magma_file_one_map(files, outFile, field)





def call_magma():




ranks_dir = args.rank_dir;
if not os.path.isdir(ranks_dir):
    os.makedirs(ranks_dir)

outd = ranks_dir + "/rank_" + (args.map_dir).split("/")[-1] 

make_magma_file_from_dir(args.map_dir, outd, args.field)
