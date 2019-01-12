#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args=sys.argv[1:]

if len(args) != 2:
	print(f"USAGE: {sys.argv[0]} <reads1> <reads2>",file=sys.stderr)
	print(f"DESCRIPTION: Takes two files of lists of reads, and merges them into a paired list file.", file=sys.stderr)
	sys.exit(1)

reads1=args[0]
reads2=args[1]

# Read both files into their individual lists
with open(reads1,"r") as f:
	r1=list(f)

with open(reads2,"r") as f:
	r2=list(f)

# Iterate through both lists simultaneously and print
for x,y in zip(r1,r2):
	line=x.rstrip() + " " + y.rstrip()
	print(line)
# Print to STDOUT, redirected to file in bash script


