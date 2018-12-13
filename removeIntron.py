#!/usr/bin/env python
from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) != 5:
	print(f"USAGE: {sys.argv[0]} <mRNA start> <intron start> <intron end> <mRNA end> <FASTA file>")
	sys.exit(1)

RNAstart=args[0]
RNAend=args[3]
INTRONstart=args[1]
INTRONend=args[2]
inputFilename=args[4]

with open(inputFilename, "r") as f:
	fasta=list(f)

header=fasta[0].split(" ")[0]
sequence=fasta[1].rstrip()

exon1=sequence[RNAstart-1:INTRONstart-1]
exon2=sequence[INTRONend:RNAend]

length=INTRONend-INTRONstart+1

print(f"{header} intron length: {length}")
print(exon1 + exon2)

