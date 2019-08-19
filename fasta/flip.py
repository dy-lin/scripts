#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args=sys.argv[1:]

inputFileName=args[0]
delimiter=args[1]

if delimiter == "commas":
	sep=","
elif delimiter == "tabs":
	sep="\t"
elif delimiter == "pipes":
	sep="|"
else:
	print("Invalid delimiter", file=sys.stderr)
	sys.exit(1)

inputMatrix=[]
with open(inputFileName,"r") as f:
	for line in f:
		inputMatrix.append(line.rstrip().split(sep))

numRows=len(inputMatrix)
numCols=len(inputMatrix[0])

outputMatrix=[]

for i in range(numCols):
	column=[]
	for j in range(numRows):
		column.append(inputMatrix[j][i])
	outputMatrix.append(column)

for i in outputMatrix:
	print(sep.join(i))





