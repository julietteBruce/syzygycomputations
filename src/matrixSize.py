import argparse
import os
import os.path
import subprocess
import re



argparser = argparse.ArgumentParser();
argparser.add_argument('matrix_file')

args = argparser.parse_args()


mat_file = args.matrix_file

with open(mat_file,"r") as mF:
    sM = mF.readlines()

entries = list(map(lambda s: list(map(lambda i: int(i), s.split())), sM ))

rows = set(map(lambda s:s[0], entries))
columns = set(map(lambda s:s[1], entries))

print(len(rows))
print("\n")
print(len(columns))
