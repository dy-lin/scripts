#!/bin/bash
set -euo pipefail 
# Takes a transcript FASTA file as input

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <FASTA file>" 1>&2
	echo "DESCRIPTION: Takes a FASTA file of transcripts, and write a transcript-to-gene CSV file." 1>&2
	exit 1
fi

fasta=$1

awk 'BEGIN{OFS=",";print "Transcript_ID,Gene_ID"}/^>/ {gene=substr($1,2);sub(/-R./, "", $1); var=substr($1,2);print gene,var}' $fasta > deseq/tx2gene.csv


# Takes a file of transcripts/proteins and creates a 'dictionary' csv file linking transcripts to genes
# echo "Transcript_ID,Gene_ID" > tx2gene.csv
# if [[ "$#" -ne 1 && -p "/dev/stdin" ]]
# then
# 	stdin=$(cat)
# 	for i in $stdin
# 	do
# 		gene=$(echo "$i" | sed 's/-R.\+$//')
# 		echo "$i,$gene" >> tx2gene.csv
# 	done
# elif [[ "$#" -eq 1 && "$1" != "/dev/fd/63" ]]
# then
# 	while read line
# 	do
# 		gene=$(echo "$line" | sed 's/-R.\+$//')
# 		echo "$line,$gene" >> tx2gene.csv
# 	done < $1
# elif [[ "$#" -eq 1 && "$1" == "/dev/fd/63" ]]
# then
# 	cat $1 > temp.txt
# 	while read line
# 	do	
# 		gene=$(echo "$line" | sed 's/-R.\+$//')
# 		echo "$line,$gene" >> tx2gene.csv
# 	done < temp.txt
# 
# 	rm temp.txt
# else
# 	echo "USAGE: $(basename $0) <FILE LIST>" 1>&2
# 	echo "DESCRIPTION: Takes a list of transcripts/proteins and writes a csv dictionary." 1>&2
# 	rm tx2gene.csv
# 	exit 1
# fi
# 

