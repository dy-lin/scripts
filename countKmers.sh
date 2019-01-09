#!/bin/bash

if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <FASTQ/A file> <K-mer size>" 1>&2
	echo "DESCRIPTION: Takes a FASTQ/A file and k-mer size and outputs the number of unique k-mers." 1>&2
	exit 1
fi

file=$1
kmer=$2
processed=temp
ext=${file##*.}
zip=${file#*.}
# If <() is used
if [[ ! -e "$file" ]]
then
	echo "File not found." 1>&2
	exit 1
fi

if [[ "$file" == "/dev/fd/63" ]]
then
	first=$(head -c 1 $file)
	if [[ "$first" == "@" ]]
	then
		ext="fq"
	elif [[ "$first" == ">" ]]
	then
		ext="fa"
	else
		echo "Unsupported file type." 1>&2
		exit 1
	fi
fi

# Find file type, as to extract sequences only
if [[ "$ext" == "fa" || "$ext" == "fasta" ]]
then
	awk '!/^>/ {print}' $file > $processed
elif [[ "$ext" == "fq"|| "$ext" == "fastq"  ]]
then
	awk '/^@/ {x=NR+1;next}(NR<=x){print}' $file > $processed
elif [[ "$zip" == "fq.gz" || "$zip" == "fastq.gz" ]]
then
	awk '/^@/ {x=NR+1;next}(NR<=x){print}' <(gunzip -c $file) > $processed
elif [[ "$zip" == "fa.gz" || "$zip" == "fasta.gz" ]] 
then
	awk '!/^>/ {print}' <(gunzip -c $file) > $processed
elif [[ "$zip" == "fa.tar.gz" || "$zip" == "fasta.tar.gz" ]]
then
	for i in $(tar -xvf $file)
	do
		awk '!/^>/ {print}' $i
	done > $processed
elif [[ "$zip" == "fq.tar.gz" || "$zip" == "fastq.tar.gz" ]]
then
	for i in $(tar -xvf $file)
	do
		awk '/^@/ {x=NR+1;next}(NR<=x){print}' $i
	done > $processed
else
	echo "Unsupported file type." 1>&2
	exit 1
fi

countKmers.py $processed $kmer

# Remove temp file created
rm $processed
