#!/bin/bash

if [ "$#" -ne 3 ]
then
	echo "USAGE: $(basename $0) <START> <END> <FASTA file>"
	exit 1
fi

if [ "$1" -gt "$2" ]
then
	slice.py $1 $2 <(seqtk seq -r $3)
else
	slice.py $1 $2 <(seqtk seq $3)
fi
