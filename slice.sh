#!/bin/bash

if [ "$#" -ne 3 ]
then
	echo "USAGE: $(basename $0) <START> <END> <FASTA file>"
	exit 1
fi

slice.py $1 $2 <(seqtk seq $3)
