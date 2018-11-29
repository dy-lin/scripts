#!/bin/bash

if [ "$#" -ne 3 ]
then
	echo "USAGE: $(basename $0) <START> <END> <FASTA file>"
	exit 1
fi

filename=$(basename $3)
extension="${filename##*.}"
filename=$(basename $3 $extension)

seqtk seq $3 > ${filename}seqtk.fa
slice.py $1 $2 ${filename}seqtk.fa
rm ${filename}seqtk.fa
