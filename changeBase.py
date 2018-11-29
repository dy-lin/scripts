#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) != 4:
	print("Usage: %s <Input File> <Base Position> <Original Base> <New Base>" % sys.argv[0])
	sys.exit(1)

directory=os.getcwd()

inputFileName=args[0]
pos=int(args[1])
oldbase=args[2].upper()
newbase=args[3].upper()
bases=["A", "C", "G", "T"]

if (oldbase not in bases) or (newbase not in bases):
	print("Invalid base.")
	sys.exit(1)

with open(inputFileName, "r") as f:
	inputFile=list(f)

header=inputFile[0]
sequence=header[1:-1]
fasta=list(inputFile[1].rstrip())

outputFileName=""

if "/dev/fd/" in inputFileName:
	outputFileName = directory + "/" + sequence + ".pos" + str(pos) + "." + oldbase + "to" + newbase + ".fa"
elif "." in inputFileName:
	filename=inputFileName.split(".")
	extLen=len(filename[-1])+1
	outputFileName = inputFileName[:-extLen] + ".pos" + str(pos) + "." + oldbase + "to" + newbase + ".fa"
else:
	outputFileName = inputFileName +  ".pos" + str(pos) + "." + oldbase + "to" + newbase + ".fa"

#Position 1 specified by user is index 0
index=pos-1

if fasta[index] != oldbase:
	print("The position specified does not match the base specified.")
	sys.exit(1)
else: 
	fasta[index] = newbase
	with open(outputFileName, "w") as f:
		f.write(header)
		for letter in fasta:
			f.write(letter)
	print(f"{oldbase} changed to {newbase} at position {pos}.")

