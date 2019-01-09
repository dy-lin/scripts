#!/usr/bin/env python
from __future__ import print_function
import os
import sys

# Given mRNA and intron coordinates (manually extracted from a GFF file), print a FASTA file to screen with intron segment removed
args = sys.argv[1:]

if len(args) != 5:
	print(f"USAGE: {sys.argv[0]} <mRNA start> <intron start> <intron end> <mRNA end> <FASTA file>",file=sys.stderr)
	print(f"DESCRIPTION: Takes mRNA and intron coordinates with the corresponding FASTA file (scaffold) and removes the introns so as to output a spliced mRNA transcript.",file=sys.stderr)
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

# Should an actual file be desired, use bash redirection to save the output as a file.

print(f"{header} intron length: {length}")
print(exon1 + exon2)

