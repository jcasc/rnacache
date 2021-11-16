#!/usr/bin/env python3

import sys
import pysam
from collections import defaultdict
import json
from multiprocessing import Pool
import os



def compare_sets(lhs, rhs):
    if len(lhs) == 0 and len(rhs) == 0:
        return "nomap"
    elif len(lhs) == 0: # len(rhs) > 0
        return "reject"
    elif len(rhs) == 0: # len(lhs) > 0
        return "uniq"

    if len(lhs) < len(rhs):
        inter = sum(x in rhs for x in lhs)
    else:
        inter = sum(x in lhs for x in rhs)

    if len(lhs) == len(rhs):
        if inter == len(lhs):
            return "eq"

    elif len(lhs) < len(rhs):
        if inter == len(lhs):
            return "sub"
    
    else: # len(lhs) > len(rhs)
        if inter == len(rhs):
            return "super"
    
    if inter == 0:
        return "disj"
    else:
        return "inter"

num_reads = int(sys.argv[1])

tid = defaultdict(lambda:len(tid))
filenames = sys.argv[3:]
for x in filenames:
    if not x.endswith(".bam"):
        print("Error: Expected .bam ending for", x, file=sys.stderr)
        sys.exit(1)

mcmap = {}

with open(sys.argv[2], "r") as infile:
    print("mcr16", file=sys.stderr)
    for line in infile:
        if line.startswith("#"):
            continue
        if len(mcmap) % 1000000 == 0:
            print(len(mcmap), file=sys.stderr)
        splt = line.split("\t")
        if len(splt[2]) > 0:
            qname = splt[0][:-2] if splt[0][-2] == "/" else splt[0]
            mcmap[qname] = set(tid[x.split(":")[0]] for x in splt[2].split(","))
            # print(qname, file=sys.stderr)

def process_file(filename):
    print(filename, file=sys.stderr)
    result = defaultdict(int)
    mapping = set()
    qname = None
    with pysam.AlignmentFile(filename, "rb") as alnFile:
        for count, rec in enumerate(alnFile):
            if count % 1000000 == 0:
                print(sum(result.values()), file=sys.stderr)
            
            qnext = rec.query_name[:-2] if rec.query_name[-2]=="/" else rec.query_name # kallisto fix
            # print("QNEXT", qnext, file=sys.stderr)
            # print("QNAME", qname, result[filename], file=sys.stderr)
            if qnext != qname:
                if len(mapping)>0:
                    mcset = mcmap[qname] if qname in mcmap else ()
                    result[compare_sets(mcset, mapping)] += 1
                    mapping.clear()
                qname = qnext
                # print("IF QNAME", qname, result[filename], file=sys.stderr)
            
            if not rec.is_unmapped:
                mapping.add(tid[rec.reference_name.split("|")[0]])

        if qname != None:
            mcset = mcmap[qname] if qname in mcmap else ()
            result[compare_sets(mcset, mapping)] += 1

    result["nomap"] = num_reads-len(mcmap)-result["reject"]
    result["uniq"] = num_reads-sum(result.values())
    return result

with Pool(len(filenames)) as p:
    results = p.map(process_file, filenames)
result = {os.path.split(k)[1].split('.bam')[0]:v for k, v in zip(filenames, results)}
    
print(json.dumps(result))
