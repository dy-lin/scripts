#!/bin/bash
set -eu -o pipefail
# Extract intron lengths given a GFF file
# Prints results to screen, use bash redirection to save file to a tsv

if [ "$#" -ne 1 ]
then
	echo "USAGE: $(basenmae $0) <GFF file>" 1>&2
	echo "DESCRIPTION: Takes a GFF file, extracts the intron features and outputs the intron lengths." 1>&2
	exit 1
fi
gff=$1

IFS=$'\n'

for line in $(cat $gff)
do
	feature=$(echo $line | awk -F "\t" '{print $3}')
	if [ "$feature" == "mRNA" ]
	then
		name=$(echo $line | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//')
	fi
	if [ "$feature" == "intron" ]
	then
		begin=$(echo $line | awk -F "\t" '{print $4}')
		end=$(echo $line | awk -F "\t" '{print $5}')
		length=$((end-begin+1))
		echo -e "$name\t$length"
	fi
done
	
