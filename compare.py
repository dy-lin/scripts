#!/usr/bin/env python
from __future__ import print_function
import os
import sys

#Take in a FASTA file and compare sequences

args = sys.argv[1:]

if len(args) != 1:
	print(f"Usage: {sys.argv[0]} <FASTA file>")
	sys.exit(1)

infile = args[0]
with open(infile,"r") as f:
	fasta = list(f)

seqID = [i.split(" ")[0][1:].rstrip() for i in fasta[0::2]]
seqs = [i.rstrip() for i in fasta[1::2]]

length = len(seqs[0])

for i in seqs:
	if len(i) != length:
		print("Your sequences are not of the same length and cannot be directly compared.")
		sys.exit(1)
	if len(seqID) != len (seqs):
		print("There is a different number of sequence IDs than sequences.")
		sys.exit(1)

listLen = len(seqs)

matched = False
for i in range(listLen-1):
	for j in range(i+1,listLen):
		if (seqs[i] == seqs[j]):
			matched = True
			print(f"{seqID[i]} is identical to {seqID[j]}.")

if matched == False:
	print("All sequences are unique!")
