#!/usr/bin/env python
from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) !=2:
	print("Usage: %s <partial sequence> <FASTA file>" % sys.argv[0])
	sys.exit(1)

substring=args[0]
filename=args[1]

with open(args[1],"r") as f:
	seq = list(f)

sequence = seq[1].rstrip()

start=sequence.index(substring)+1
end=len(substring) + start-1

print(f"{start} to {end}")



