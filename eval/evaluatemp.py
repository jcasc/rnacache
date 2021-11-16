#!/usr/bin/env python3

import pysam
import sys
from multiprocessing import Pool
import os.path

num_reads = int(sys.argv[1])

def get_refname(alnFile, r):
    return alnFile.getrname(r).split('|',1)[0]

def get_origin(q):
    splt = q.split(':',3)
    return splt[2] if len(splt)>=3 else ""

def evaluate(filename):
    with pysam.AlignmentFile(filename, 'r') as alnFile:
        total_matches = 0
        origin_mapped = 0
        reads_aligned = 0
        origin_weighted = 0

        query = ""
        origin = ""
        matches = set()
        for rec in alnFile:
            query_name = rec.qname[:-2] if rec.qname[-2]=="/" else rec.qname # kallisto fix
            if query_name != query:
                if len(matches)>0:
                    reads_aligned += 1
                    total_matches += len(matches)
                    if origin in matches:
                        origin_mapped += 1
                        origin_weighted += 1/len(matches)
                    matches.clear()
                query = query_name
                origin = get_origin(query)
            if not rec.is_unmapped:
                rname = get_refname(alnFile, rec.rname)
                matches.add(rname)
        
        if len(matches)>0:
            reads_aligned += 1
            total_matches += len(matches)
            if origin in matches:
                origin_mapped += 1
                origin_weighted += 1/len(matches)


    with open(f"{prefix}/{os.path.split(filename)[1][:-4]}_eval", "w") as outfile:
        # print(reads_aligned, origin_mapped, origin_weighted, total_matches, file=outfile)
        # print(reads_aligned/num_reads, origin_mapped/num_reads, origin_weighted/reads_aligned, total_matches/num_reads, file=outfile)

        print(f"# Total reads (or pairs): {num_reads}\n"
              f"# Total matches:          {total_matches}\n"
              f"# Origin found:           {origin_mapped}\n"
              f"# Reads aligned:          {reads_aligned}\n"
              f"# Recall:                 {origin_mapped/num_reads}\n"
              f"# Hits per read:          {total_matches/num_reads}\n"
              f"# Hits per aligned read:  {total_matches/reads_aligned}\n"
              f"# Total True Hit Rate:    {origin_mapped/total_matches}\n"
              f"# Mean True Hit Rate:     {origin_weighted/reads_aligned}\n", file=outfile)

prefix = sys.argv[2]
filenames = sys.argv[3:]
if any(not x.endswith(".bam") for x in filenames):
    print("files must have .bam ending!", file=sys.stderr)
else:
    with Pool(len(filenames)) as p:
        p.map(evaluate, filenames)
