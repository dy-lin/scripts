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
header = seq[0][1:].split(" ")[0].rstrip()
sequence = seq[1].rstrip()
count = sequence.count(substring)

i = 1
pointer=0
tempSeq = sequence
print(f"Searching for {substring} in {header}...")
print(f"Number of matches: {count}")
start=0
end=0
try:
	while i <= count:
		if i==1:
			start=tempSeq.index(substring)
			end=len(substring) + start-1
			print(f"{start+1} to {end+1}")
		else:
			start=tempSeq.index(substring) + end +1
			end=len(substring) + start-1
			print(f"{start+1} to {end+1}")
		i+=1
		tempSeq=sequence[end+1:]
except ValueError as error:
	print("Partial sequence not found.")
