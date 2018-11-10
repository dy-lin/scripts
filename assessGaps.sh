#!/bin/bash
if [ "$#" -eq 0 ]
then
	echo "USAGE: $(basename $0) <FASTA file> [GFF file]"
	exit 1
fi

filename=$(basename $1)
extension="${filename##*.}"
filename=$(basename $1 $extension)

seqtk seq $1 > ${filename}seqtk.fa
if [ "$#" -eq 2 ]
then
	findGaps.py ${filename}seqtk.fa $2
else
	findGaps.py ${filename}seqtk.fa
fi

if [ -e "${filename}.seqtk.gaps.tsv" ]
then
	column -t -s$'\t' ${filename}seqtk.gaps.tsv
fi
