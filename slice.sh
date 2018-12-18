#!/bin/bash

if [ "$#" -ne 3 ]
then
	echo "USAGE: $(basename $0) <START> <END> <FASTA file> or"
	echo "USAGE: $(basename $0) <POSITION> <STRAND> <FASTA file>, where STRAND is + or -"
	exit 1

fi
if [ "$2" == "-" ] || [ "$2" == "+" ]
then
	slice.py $1 $2 $3
elif [ "$1" -gt "$2" ]
then
	slice.py $1 $2 <(seqtk seq -r $3)
elif [ "$1" -lt "$2" ]
then
	slice.py $1 $2 <(seqtk seq $3)
else
	echo "START and END cannot be the same. If one base at a specific position is needed:"
	echo "USAGE: $(basename $0) <POSITION> <STRAND> <FASTA file>, where STRAND is + or -"
fi

