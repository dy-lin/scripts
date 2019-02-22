#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args = sys.argv[1:]

if len(args) != 2 and len(args) !=1:
	print(f"USAGE: {sys.argv[0]} <FASTA file> [GFF file]",file=sys.stderr)
	print(f"DESCRIPTION: Takes a FASTA file and outputs its gaps and a scaftig file. If a GFF file is given, features containing gaps will also be listed.",file=sys.stderr)
	sys.exit(1)

file = args[0]

with open(file, "r") as f:
	fileList=list(f)

# Isolate header with just gene name, no trailing information, such as eAED scores, etc.
header=fileList.pop(0).split(" ")[0]
sequence=fileList[0].rstrip()
gap = False
start = 0

# Stores start and end STRINGS of gaps
gaps = []

# Stores length of each GAP
lengths=[]

# Stores coordinates of each gap (Start and end as integers)
borders=[]
for index, base in enumerate(sequence):
	if gap == False:
		if base.upper() == "N":
			gap = True
			start = index+1
	else:
		if base.upper() !="N":
			gap = False
			end = index
			if start == end:
				 gaps.append(str(start))
			else:
				gaps.append("" + str(start) + " to " + str(end))
			borders.append([start,end])
			lengths.append(end-start+1)

# The largest gap length
if len(lengths) != 0:
	maxgap=max(lengths)
else:
	maxgap=0
outFile=[]
num_gaps = len(gaps)

# Remove leading directory string, keeping only filename
filename=file.split("/")[-1]

# Formatting Gap report headings and columns
outFile.append("Gaps in " + filename + "\n")
outFile.append("Number of gaps: " + str(num_gaps) + "\n")
outFile.append("Longest gap length: " + str(maxgap) + "\n\n")
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
	# If a GFF file was given, show which features contain gaps
	if len(args) == 2:
		gff = args[1]
		with open(gff, "r") as f:
			gffFile=list(f)
		
		# Storing a list of lists
		# Storing a list of a list of each field in the GFF
		gffList=[]
		for index, line in enumerate(gffFile):
			gffList.append(line.rstrip().split("\t"))
		# Storing columns/fields of 3, 4, 5, and parts of 8
		features=[]
		for i in gffList:
			if i[0][0] != "#":
				temp=i[8].split(";")
				features.append([i[2],i[3],i[4],temp[0][3:]])
	
		# Storing features that are not mRNA or genes, as to avoid repeat 'features' with gaps in the report
		nonmRNA=[]
		for i in features:
			if i[0] != "mRNA" and i[0] != "gene":
				nonmRNA.append(i)
		
		# Storing gene features that contain gaps
		genes=[]	
		for index, pair in enumerate(borders):
			for i in range(pair[0],pair[1]+1):
				for j in nonmRNA:
					if i in range(int(j[1]),int(j[2])+1):
						genes.append([gaps[index], j[3]])

		# Populating scaftigs list with one element - the first scaftig from 1 to border of first gap-1
		# Then populating the rest using the same method, but a special case for the last scaftig as with the first
		scaftigs=[[1,borders[0][0]-1]]
		for i in range(len(borders)):
			if i == len(borders)-1:
				scaftigs.append([borders[i][1]+1, len(sequence)])
			else:
				scaftigs.append([borders[i][1]+1,borders[i+1][0]-1])

		# Storing which gene features are on which scaftigs
		feats=[]
		for index, pair in enumerate(scaftigs):
			for i in range(pair[0],pair[1]+1):
				for j in nonmRNA:
					if i in range(int(j[1]),int(j[2])+1):
						if [str(pair[0]) + " to " + str(pair[1]),str(j[3])] not in feats:
							feats.append([str(pair[0]) + " to " + str(pair[1]), str(j[3])])

		with open(outFileName, "a") as f:
			f.write("\nGFF: " + gff + "\n")
			f.write("Gap Position\tGene:Feature\n")
			for i in genes:
				f.write(i[0]+ "\t" + i[1] + "\n")
			f.write("\nScaftig Position\tGene:Feature\n")
			for i in feats:
				f.write(i[0] + "\t" + i[1] + "\n")

	# If there is no GFF given, skip the gene feature parts
	else:
		scaftigs=[[1,borders[0][0]-1]]
		for i in range(len(borders)):
			if i == len(borders)-1:
				scaftigs.append([borders[i][1]+1, len(sequence)])
			else:
				scaftigs.append([borders[i][1]+1,borders[i+1][0]-1])
		with open(scaftigsOUT, "w") as f:
			for seq in scaftigs:
				f.write(header.rstrip() + "-" + str(seq[0]) + ":" + str(seq[1]) + " Length: " + str(seq[1]-seq[0]+1) +  "\n")
				slice = sequence[seq[0]-1:seq[1]]
				f.write(slice + "\n")

	# Mini report printed to screen
	print("Gaps in " + filename)
	print("Number of gaps: " + str(num_gaps))
	print("Longest gap length: " + str(maxgap))
else:
	print("Gaps in " + filename)
	print("Number of gaps: " + str(num_gaps))
