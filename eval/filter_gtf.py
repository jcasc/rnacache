#!/usr/bin/env python3

import sys

transcript_file = sys.argv[1]
in_file = sys.argv[2]

with open(transcript_file, "r") as f:
    ids = set(line.split("|", 1)[0][1:] for line in f if line.startswith(">"))

with open(in_file, "r") as f:
    for line in f:
        if line[0]=='#':
            continue
        tmp = line.split()
        if tmp[10]!="transcript_id" or tmp[11][1:-2] in ids:
            print(line,end="")
