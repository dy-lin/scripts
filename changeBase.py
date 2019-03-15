#!/usr/bin/env python

from __future__ import print_function
import os
import sys

# Used when Pilon misses a base that is not in consensus with the read alignments

args = sys.argv[1:]

if len(args) != 4:
	print(f"USAGE: {sys.argv[0]} <Input File> <Base Position> <Original Base> <New Base>",file=sys.stderr)
	print(f"DESCRIPTION: Takes a FASTA file, a base position, the original base and the new base as input, outputting the newly modified FASTA file.",file=sys.stderr)
	sys.exit(1)

# Get current working directory
directory=os.getcwd()

inputFileName=args[0]

# Position of the change
pos=int(args[1])

# Original Base to be changed
oldbase=args[2].upper()

# New base; base to be
newbase=args[3].upper()

# List of acceptable bases
bases=["A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]

# Print error message if any of the arguments are invalid bases
if (oldbase not in bases) or (newbase not in bases):
	print("Invalid base.")
	sys.exit(1)

with open(inputFileName, "r") as f:
	inputFile=list(f)

header=inputFile[0]

# Sequence is header without the trailing newline and the first >
sequence=header[1:-1]

# The FASTA sequence without the trailing newline
fasta=list(inputFile[1].rstrip())

outputFileName=""

# Naming the output file to reflect positional base changes in the name
if "/dev/fd/" in inputFileName:
	outputFileName = directory + "/" + sequence + ".pos" + str(pos) + "." + oldbase + "to" + newbase + ".fa"
elif "." in inputFileName:
	filename=inputFileName.split(".")
	extLen=len(filename[-1])+1
	outputFileName = inputFileName[:-extLen] + ".pos" + str(pos) + "." + oldbase + "to" + newbase + ".fa"
else:
	outputFileName = inputFileName +  ".pos" + str(pos) + "." + oldbase + "to" + newbase + ".fa"

# Position 1 specified by user is index of 0
index=pos-1

# Error message if the position given and the original base do not corresond in the file
# If it corresponds, change the base to the specified new base
if fasta[index] != oldbase:
	print("The position specified does not match the base specified.")
	sys.exit(1)
else: 
	fasta[index] = newbase
	# Write new sequence out into the new output FASTA file.
	with open(outputFileName, "w") as f:
		f.write(header)
		for letter in fasta:
			f.write(letter)
	print(f"{oldbase} changed to {newbase} at position {pos}.")
