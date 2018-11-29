#!/usr/bin/env python
from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) != 2 and len(args) !=1:
	print("Usage: %s <FASTA file> [GFF file]" % sys.argv[0])
	sys.exit(1)

file = args[0]
#file = input("Filename: ")
with open(file, "r") as f:
	fileList=list(f)

header=fileList.pop(0).split(" ")[0]
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
if num_gaps !=0:
	outFile.append("Gap Number\tGap Position\tGap Length\n")
	for i,N in enumerate(gaps):
		outFile.append(str(i+1) + "\t" + N + "\t" + str(lengths[i])+"\n")
	scaftigsOUT=""
	outFileName=""
	if "." in file:
		name=file.split(".")
		extLen=len(name[-1])+1
		outFileName=file[:-extLen] + ".gaps.tsv"
		scaftigsOUT=file[:-extLen] + ".scaftigs.fa"
	else:
		outFileName=file + ".gaps.tsv"
		scaftigsOUT=file + ".scaftigs.fa"
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
			if i[0] != "mRNA" and i[0] != "gene":
				nonmRNA.append(i)
		genes=[]
		
		for index, pair in enumerate(borders):
			for i in range(pair[0],pair[1]+1):
				for j in nonmRNA:
					if i in range(int(j[1]),int(j[2])+1):
						genes.append([gaps[index], j[3]])
		
		scaftigs=[[1,borders[0][0]-1]]
		for i in range(len(borders)):
			if i == len(borders)-1:
				scaftigs.append([borders[i][1]+1, len(sequence)])
			else:
				scaftigs.append([borders[i][1]+1,borders[i+1][0]-1])
		
		feats=[]
		for index, pair in enumerate(scaftigs):
			for i in range(pair[0],pair[1]+1):
				for j in nonmRNA:
					if i in range(int(j[1]),int(j[2])+1):
						if [str(pair[0]) + " to " + str(pair[1]),str(j[3])] not in feats:
							feats.append([str(pair[0]) + " to " + str(pair[1]), str(j[3])])
	#	feats_unique=set(feats)
	#	feats_unique=sorted(feats_unique)
		with open(outFileName, "a") as f:
			f.write("\nGFF: " + gff + "\n")
			f.write("Gap Position\tGene:Feature\n")
			for i in genes:
				f.write(i[0]+ "\t" + i[1] + "\n")
			f.write("\nScaftig Position\tGene:Feature\n")
			for i in feats:
				f.write(i[0] + "\t" + i[1] + "\n")

	else:
		scaftigs=[[1,borders[0][0]-1]]
		for i in range(len(borders)):
			if i == len(borders)-1:
				scaftigs.append([borders[i][1]+1, len(sequence)])
			else:
				scaftigs.append([borders[i][1]+1,borders[i+1][0]-1])
		with open(scaftigsOUT, "w") as f:
			for seq in scaftigs:
				f.write(header + "-" + str(seq[0]) + ":" + str(seq[1]) + " Length: " + str(seq[1]-seq[0]+1) +  "\n")
				slice = sequence[seq[0]-1:seq[1]]
				f.write(slice + "\n")

else:
	print("Gaps in " + file)
	print("Number of gaps: " + str(num_gaps))
