#!/bin/bash
stdin=$(cat)
PROGRAM=$(basename $0)
if [[ ! -p "/dev/stdin" && "$#" -lt 1 ]]
then
	echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) as input, and outputs each sequence in a separate FASTA file." 1>&2
	exit 1
fi

# Read ID and sequence
if [[ "$#" -ge 1 ]]
then
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
elif [[ -p "/dev/stdin" ]]
then
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
	done < <(echo "$stdin")
fi

