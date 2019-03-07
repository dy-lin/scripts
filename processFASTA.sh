#!/bin/bash

PROGRAM=$(basename $0)

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Processes a multi-line FASTA file into one line." 1>&2
	exit 1
fi

num=1

for fasta in $@
do

	if [[ "$fasta" != *.fa && "$fasta" != *.fasta ]]
	then
		echo "ERROR: Input file is not a FASTA file." 1>&2
		echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
		echo "DESCRIPTION: Processes a multi-line FASTA file into one line." 1>&2
		exit 1
	fi

	if [[ -e "$fasta" ]]
	then
		echo "Processing $(basename $fasta)..." 1>&2
	else
		echo "$(basename $fasta) does not exist!" 1>&2
		continue
	fi

	path=$(dirname $fasta)
	# In case there are other scripts using temp files of similar names.
	if [[ -e "${path}/temp${num}.fa" ]]
	then
		num=$((num+1))
	fi
	seqtk seq $fasta > ${path}/temp${num}.fa
	mv ${path}/temp${num}.fa $fasta
done

