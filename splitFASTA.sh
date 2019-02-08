#!/bin/bash
PROGRAM=$(basename $0)

if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) as input, and outputs each sequence in a separate FASTA file." 1>&2
	exit 1
fi

# Read ID and sequence
while read id
do
	filename=$(echo $id | awk '{print $1}' | sed 's/^>//').fa
	read seq
	# If output file already exists, there are identical sequences in the file(s) - only write one.
	if [[ ! -e "$filename" ]]
	then
		echo $id >> $filename
		echo $seq >> $filename
	fi
done < <(cat $*)
