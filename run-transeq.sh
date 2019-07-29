#!/bin/bash

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <nucleotide sequence(s)>" 1>&2
	exit 1
fi

count=1
filenum=1
outfile=temp.fa
while [[ -e "$outfile" ]]
do
	filenum=$((filenum+1))
	outfile="temp${filenum}.fa"
done

for seq in $@
do
	echo -e ">header$count\n$seq" >> $outfile
	count=$((count+1))
done

transeq -sequence $outfile -outseq /dev/stdout 2> /dev/null | sed 's/X$//' | seqtk seq -
rm $outfile

