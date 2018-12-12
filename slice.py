#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) != 3:
	print("Usage: %s <Start> <End> <Input File>" % sys.argv[0])
	sys.exit(1)

x=int(args[0])
y=int(args[1])
#BLAST has one-start counting vs. slicing the string, which has 0-based counting
reverse = False
if x > y:
	# if start coordinate is greater than end coordinate, than switch them
	x,y=y,x
	reverse = True
if x == y:
	print("Start and End coordinates are the same!")
	sys.exit(1)
start=x-1
end=y-1
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


#with open(outputFileName, "w") as f:
for line in inputFile:
	if line[0] == ">":
	#	f.write(line.strip() + " Start: " + str(start+1) + " End: " + str(end+1) + "\n")
		print(line.strip() + " Start: " + str(start+1) + " End: " + str(end+1))

	else:
		#if coordinates were flipped, the reverse complement of the strand was taken, thereby changing the coordinate placement
		# to solve this, print the sequence backwards so the coordinates are correct but the bases are still in reverse order
		# then print the sliced portion in reverse, leaving the 5' -> 3' correct sliced portion with correct bases
		if reverse == True:
			rc=line.rstrip()[::-1]
			slice=rc[start:end+1]
			print(slice[::-1])
		elif slice != "":
			slice = line[start:end+1].rstrip()
	#		f.write(slice + "\n")
			print(slice)
		else:
			print("Incompatible coordinates.")
