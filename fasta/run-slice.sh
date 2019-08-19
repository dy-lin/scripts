#!/bin/bash
# set -eu -o pipefail
# Slices the FASTA at the given coordinates
# If given one coordinate, strand must be specified

if [[ "$#" -ne 3 && "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <START> <END> <FASTA file> or" 1>&2
	echo "USAGE: $(basename $0) <POSITION> <STRAND> <FASTA file>, where STRAND is + or -" 1>&2
	echo "DESCRIPTION: Takes start and end coordinates and extracts that segment from a FASTA file, or takes a position and strand orientation and extracts that one base from that FASTA file." 1>&2
	exit 1

fi

if [[ "$#" -eq 2 ]]
then
	echo "Only start coordinate detected. Using sequence length as end coordinate. Slicing..." 1>&2
	length=$(seqtk comp $2 | awk '{print $2}')
	slice.py $1 $length <(seqtk seq $2)
elif [[ "$2" == "-"  ||  "$2" == "+" ]]
then
	echo "Position and strand detected. Slicing..." 1>&2
	slice.py $1 $2 <(seqtk seq $3)
elif [ "$1" -gt "$2" ]
then
	echo "Start coordinate is greater than end coordinate. Taking the reverse complement. Slicing..." 1>&2
	slice.py $1 $2 <(seqtk seq -r $3)
elif [ "$1" -lt "$2" ]
then
	echo "Start and end coordinates detected. Slicing..." 1>&2
	slice.py $1 $2 <(seqtk seq $3)
else
	echo "START and END cannot be the same. If one base at a specific position is needed:" 1>&2
	echo "USAGE: $(basename $0) <POSITION> <STRAND> <FASTA file>, where STRAND is + or -" 1>&2
fi

