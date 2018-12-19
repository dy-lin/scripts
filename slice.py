#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args = sys.argv[1:]
#This part unnecessary as it is in the bash script that calls this python script
#if len(args) != 3:
#	print("Usage: %s <Start> <End> <Input File>" % sys.argv[0])
#	sys.exit(1)

x=int(args[0])
y=args[1]
inputFileName=args[2]

#BLAST has one-start counting vs. slicing the string, which has 0-based counting
neg = False
same = False
reverse = False

#Script is compatible if only one coordinate is given
if y == "+":
	same = True
	neg = False
	print("One position detected. Returning the base at that positon.")
	base=x-1
elif y == "-":
	same = True
	neg = True
	print("One position detected. Returning the base at that position.")
	base=x-1
elif x > int(y):
	# if start coordinate is greater than end coordinate, than switch them
	reverse = True
	start=int(y)-1
	end=x-1
else:
	start=x-1
	end=int(y)-1

with open(inputFileName, "r") as f:
	inputFile=list(f)
#Uncomment to write to file in python script
#Easier to use bash redirection to save file 

#outputFileName=""
#Find extension

#length=end-start+1
#if "." in inputFileName:
#	filename=inputFileName.split(".")
#	extLen=len(filename[-1])+1
#	outputFileName = inputFileName[:-extLen] + "." + str(start+1) + "to" + str(end+1) + "bp.fa"
#else:
#	outputFileName = inputFileName + "." + str(start+1) +  "to" +  str(end+1) + "bp.fa"


#with open(outputFileName, "w") as f:
for line in inputFile:
	if line[0] == ">":
	#	f.write(line.strip() + " Start: " + str(start+1) + " End: " + str(end+1) + "\n")
		if same == True:
			print(line.strip() + " Position: " + str(base+1))
		else:
			print(line.strip() + " Start: " + str(start+1) + " End: " + str(end+1))

	else:
		#if coordinates were flipped, the reverse complement of the strand was taken in the bash script, thereby changing the coordinate placement
		# to solve this, print the sequence backwards so the coordinates are correct but the bases are still in reverse order
		# then print the sliced portion in reverse, leaving the 5' -> 3' correct sliced portion with correct bases
	
		# if only one base is needed, easier to do give complementary base manually
		if same == True and neg == True:
			if line[base].upper() == 'A':
				print('T')
			elif line[base].upper() == 'T':
				print('A')
			elif line[base].upper() == 'C':
				print('G')
			elif line[base].upper() == 'G':
				print('C')
			else:
				print("Invalid base.")
		elif same == True and neg == False:
			print(line[base])
		elif reverse == True:
			rc=line.rstrip()[::-1]
			slice=rc[start:end+1]
			print(slice[::-1])
	# if the slice is NOT an empty string (in the case that it is if only 1 coordinate is given), follow normal procedure to slice
		elif slice != "":
			slice = line[start:end+1].rstrip()
	#		f.write(slice + "\n")
			print(slice)
		else:
			print("Incompatible coordinates.")
