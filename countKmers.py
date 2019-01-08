#!/usr/bin/env python

from __future__ import print_function
import os
import sys

args = sys.argv[1:]

filename=args[0]
klen=int(args[1])

with open(filename,"r") as f:
	raw=list(f)

# Strip newline characters
processed=[]
for seq in raw:
	processed.append(seq.rstrip())

kmers=[]

# Split into K-mers, only storing those of valid length
for seq in processed:
	for i in range(klen):
		for j in range(i,len(seq),klen):
			if j <= len(seq)-klen:
				kmers.append(seq[j:j+klen])

#print(f"All K-Mers: {kmers}")
#print(f"Number: {len(kmers)}")

# Turn list into a set of unique K-mers
unique=set(kmers)

#print(f"Unique K-Mers: {unique}")
#print(f"Number: {len(unique)}")
print(f"Number of Unique K-mers: {len(unique)}")
