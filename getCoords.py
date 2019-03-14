#!/usr/bin/env python
from __future__ import print_function
import os
import sys

# Given a partial sequence, print to screen the coordinates of the FASTA
# Mostly used in parsing through gmap -A results where numbering is abysmal, as well as in parsing multiple sequence alignments

args = sys.argv[1:]

if len(args) !=2:
	print(f"USAGE: {sys.argv[0]} <partial sequence> <FASTA file>",file=sys.stderr)
	print(f"DESCRIPTION: Takes a partial sequence string and aligns it to the given FASTA file, outputting the start and end coordinates of the string sequence.", file=sys.stderr)
	sys.exit(1)

substring=args[0]
filename=args[1]

with open(args[1],"r") as f:
	seq = list(f)

header = seq[0][1:].split(" ")[0].rstrip()
sequence = seq[1].rstrip().upper()
count = sequence.count(substring)

i = 1
pointer=0

# Store sequence as to not lose the original argument
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
		# Update tempSeq as to exclude the 'parsed' part of the string since 's.index(i)' only returns the first hit and not others after
		# This method will make sure all hits are printed
		# However all indicies are still kept relevant to the original string
		tempSeq=sequence[end+1:]
except ValueError as error:
	print("Partial sequence not found.") 
