#!/usr/bin/python3

import sys
import re

def mode(lst):
    lst = list(lst);
    lst.sort();
    best = None;
    bestCount = 0;
    curr = None;
    currCount = 0;
    for i in range(0,len(lst)):
        if lst[i] == curr:
            currCount+=1;
        else:
            if best==None or currCount>bestCount:
                best = curr;
                bestCount = currCount;
            curr = lst[i];
            currCount = 1;
    if best==None or currCount>bestCount:
        best = curr;
        bestCount = currCount;
    return best;

for line in sys.stdin:
    lineParts = line.split();
    md = tuple(map(int,re.findall("\d+",lineParts[0])));
    print("{} {}".format(" ".join(map(str,md)),mode(map(int,lineParts[1:]))));
