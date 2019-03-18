#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args=sys.argv[1:]

if len(args) != 1:
	print(f"USAGE: {sys.argv[0]} <FASTA file>", file=sys.stderr)
	print(f"DESCRIPTION: Takes a FASTA file, and prints 'nucl' or 'prot'.", file=sys.stderr)
	sys.exit(1)

sequence=args[0]
bases=['A', 'C', 'G', 'T', 'U' ]
aas=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
amb=['R','Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N' ]
nucl=True
aa=False
ambiguous=False
for letter in sequence:
	if letter not in bases:
		nucl=False
		if letter not in aas:
			aa=False
			if letter in amb:
				ambiguous=True
				break
			break
		aa=True
		break


if nucl == True:
	print("nucl")
if nucl == False and aa == True:
	print("prot")
if nucl == False and aa == False and ambiguous == False:
	print("invalid")
if nucl == False and aa == False and ambiguous == True:
	print("ambiguous")


	



