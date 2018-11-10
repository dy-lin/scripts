#!/usr/bin/python

from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) != 3:
	print("Usage: %s <Start> <End> <Input File>" % sys.argv[0])
	sys.exit(1)

#BLAST has one-start counting vs. slicing the string, which has 0-based counting
start=int(args[0])-1
end=int(args[1])-1
inputFileName=args[2]

with open(inputFileName, "r") as f:
	inputFile=list(f)

outputFileName=""
#Find extension

#length=end-start+1
if "." in inputFileName:
	filename=inputFileName.split(".")
	extLen=len(filename[-1])+1
	outputFileName = inputFileName[:-extLen] + "." + str(start+1) + "to" + str(end+1) + "bp.fa"
else:
	outputFileName = inputFileName + "." + str(start+1) +  "to" +  str(end+1) + "bp.fa"


with open(outputFileName, "w") as f:
	for line in inputFile:
		if line[0] == ">":
			f.write(line.strip() + " Start: " + str(start+1) + " End: " + str(end+1) + "\n")
			print(line.strip() + " Start: " + str(start+1) + " End: " + str(end+1))
		else:
			slice = line[start:end+1]
			f.write(slice + "\n")
			print(slice)
