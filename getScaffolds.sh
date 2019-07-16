#!/bin/bash

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <FASTA file> <scaffold name(s)>" 1>&2
	echo "DESCRIPTION: Gets the sequence of each given scaffold into separate FASTA files." 1>&2
	exit 1
fi

fasta=$1
shift
scaffolds=$@

for i in $scaffolds
do
	seqtk subseq $fasta <(echo "$i") > ${i}.scaffold.fa
done


