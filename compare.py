#!/usr/bin/env python
from __future__ import print_function
import os
import sys

# Take in a FASTA file and compare sequences
# Quick comparison of sequences of equal length that appear similar, on the go
# Tells you if any sequences in a given FASTA file are identical or NOT
# For better alignments, do a pairwise alignment or multiple sequence alignment
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
	if len(seqID) != len (seqs):
		print("There is a different number of sequence IDs than sequences.")
		sys.exit(1)

listLen = len(seqs)

# Sequences are only considered identical if they are equal in length and residues
# If one sequence is a partial sequence of another, they are not counted as identical because their lengths differ
matched = False
for i in range(listLen-1):
	for j in range(i+1,listLen):
	#	if (len(seqs[i]) != len(seqs[j])):
		#	print(f"{seqID[i]} is a different length than {seqID[j]} and cannot be directly compared.")
		if (seqs[i].upper() == seqs[j].upper()):
			matched = True
			print(f"{seqID[i]} is identical to {seqID[j]}.")

if matched == False:
	print("All sequences are unique!")
