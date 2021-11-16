#!/usr/bin/env python3

import sys
from math import sqrt
import json
import os
from typing import DefaultDict


def avg(a):
    return sum(a)/len(a)

def sd(a):
    if len(a) == 1:
        return 0
    mean = avg(a)
    return sqrt(sum((i-mean)**2 for i in a)/(len(a)-1))

def time(filenames, decimals=2):
    # print(filenames, sep='\n')
    times = [int(splt[0])*60 + float(splt[1]) for splt in (line.split()[-1][:-1].split('m') for f in filenames for line in open(f) if "real" in line)]
    # print(times)
    # print(avg(times), sd(times))
    # print(avg(times))
    return avg(times)

def results(filenames):
    # print(filenames, sep='\n')
    paligned, recall, hpr, hpar, tthr, mthr = [], [], [], [], [], []
    for f in filenames:
        for line in open(f):
            if "Total reads" in line:
                total = int(line.split()[-1])
            elif "Reads aligned" in line:
                aligned = int(line.split()[-1])
            elif "Hits per read" in line:
                hpr.append(float(line.split()[-1]))
            elif "Hits per aligned read" in line:
                hpar.append(float(line.split()[-1]))
            elif "Total True Hit Rate" in line:
                tthr.append(100*float(line.split()[-1]))
            elif "Mean True Hit Rate" in line:
                mthr.append(100*float(line.split()[-1]))
            elif "Recall" in line:
                recall.append(100*float(line.split()[-1]))
        paligned.append(100*aligned/total)

    longfmt="{:13.10f}"
    shortfmt="{:.1f}"
    # print(f"palgn:\t{longfmt.format(avg(paligned))}\t{shortfmt.format(avg(paligned))}\n"
    #       f"recall:\t{longfmt.format(avg(recall))}\t{shortfmt.format(avg(recall))}\n"
    #       f"hpr:\t{longfmt.format(avg(hpr))}\t{shortfmt.format(avg(hpr))}\n"
    #       f"hpar:\t{longfmt.format(avg(hpar))}\t{shortfmt.format(avg(hpar))}\n"
    #       f"tthr:\t{longfmt.format(avg(tthr))}\t{shortfmt.format(avg(tthr))}\n"
    #       f"mthr:\t{longfmt.format(avg(mthr))}\t{shortfmt.format(avg(mthr))}")
    
    return {"palgn":avg(paligned),
            "recall":avg(recall),
            "hpr":avg(hpr),
            "hpar":avg(hpar),
            "tthr":avg(tthr),
            "mthr":avg(mthr)}
    

# main


if sys.argv[1] == "time":
    command = time
    suffix = "_time"
elif sys.argv[1] == "eval":
    command = results
    suffix = "_eval"

outpath = sys.argv[2]
if os.path.exists(outpath):
    print(outpath,"already exists!", file=sys.stderr)
    sys.exit(1)

paths = sys.argv[3:]

for x in paths:
    if not x.endswith(suffix):
        print(f"ERROR: expected suffix {suffix} for {x}", file=sys.stderr)
        sys.exit(1)

groups = DefaultDict(list)
for p in paths:
    head, tail = os.path.split(p)
    groups[tail].append(head)

combos = DefaultDict(list)
for g in groups:
    combos[tuple(groups[g])].append(g)

for c in combos:
    print(c, *combos[c], len(combos[c]), "\n", sep='\n')
    

res = {g.split(suffix)[0]:command([p+"/"+g for p in groups[g]]) for g in groups}
json.dump(res, open(outpath, "w"), indent=True)

