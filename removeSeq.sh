#!/bin/bash

if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <sequence header> <FASTA file>" 1>&2
	echo "DESCRIPTION: Given a (partial) sequence header, removes those sequences from the FASTA file." 1>&2
	exit 1
fi

header=$1
file=$2
filename=${file%.*}
# For this to work, FASTA file must be in correct file format, no multi-line FASTAs.
# Remove empty lines from the original file.

sed -i '/^$/d' $file

num_seqs=$(grep -c '^>' $file)
num_lines=$(wc -l $file | awk '{print $1}')

num_correct_lines=$((num_seqs*2))

if [[ "$num_lines" -ne "$num_correct_lines" ]]
then
	seqtk seq $file > ${filename}.seqtk.fa
	file=${filename}.seqtk.fa
fi
sed -i "/$header/,+1 d" $file
