#!/usr/bin/python

import os
import sys

args = sys.argv[1:]

if len(args) != 2 and len(args) !=1:
	print("Usage: %s <FASTA file> [GFF file]" % sys.argv[0])

file = args[0]
#file = input("Filename: ")
with open(file, "r") as f:
	fileList=list(f)

fileList.pop(0)
sequence=fileList[0].rstrip()
gap = False
start = 0
gaps = []
lengths=[]
borders=[]

for index, base in enumerate(sequence):
	if gap == False:
		if base == "N":
			gap = True
			start = index+1
	else:
		if base !="N":
			gap = False
			end = index
			if start == end:
				 gaps.append(str(start))
			else:
				gaps.append("" + str(start) + " to " + str(end))
			borders.append([start,end])
			lengths.append(end-start+1)

outFile=[]
num_gaps = len(gaps)
outFile.append("Gaps in " + file + "\n")
outFile.append("Number of gaps: " + str(num_gaps) + "\n\n")

outFile.append("Gap Number\tGap Position\tGap Length\n")
for i,N in enumerate(gaps):
	outFile.append(str(i+1) + "\t" + N + "\t" + str(lengths[i])+"\n")
outFileName=file[:-2] + "gaps.tsv"
outFileString= "".join(outFile)
with open(outFileName, "w") as f:
	f.write(outFileString)

if len(args) == 2:
	gff = args[1]
	with open(gff, "r") as f:
		gffFile=list(f)
	gffList=[]
	for index, line in enumerate(gffFile):
		gffList.append(line.rstrip().split("\t"))
	features=[]
	for i in gffList:
		temp= i[8].split(";")
		features.append([i[2],i[3],i[4],temp[0][3:]])
	nonmRNA=[]
	for i in features:
		if i[0] != "mRNA":
			nonmRNA.append(i)
	genes=[]
	for index, pair in enumerate(borders):
		for i in range(pair[0],pair[1]+1):
			for j in nonmRNA:
				if i in range(int(j[1]),int(j[2])+1):
					genes.append([gaps[index], j[3]])
	with open(outFileName, "a") as f:
		
		f.write("\nGFF: " + gff + "\n")
		f.write("Gap Position\tGene:Feature\n")
		for i in genes:
			f.write(i[0]+ "\t" + i[1])
