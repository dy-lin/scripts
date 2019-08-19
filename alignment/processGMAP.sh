#!/bin/bash

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <GMAP GFF output>" 1>&2
	echo "DESCRIPTION: Takes a GMAP GFF file and outputs a tabular BLAST output file." 1>&2
	exit 1 
fi
gff=$1
#long file name, rids the .gff portion.
filename=${gff%.*}
outfile=${filename}.tsv

# print headers first
echo -e "Query\tSubject\tPID\tMatches\tMismatches\tIndels\tUnknowns\tCov\tSstart\tSend" > $outfile

#read GFF file

# take whole file without #s

IFS=$'\n'
for alignment in $(grep -v '^#' $gff | awk '/\tmRNA\t/')
do
	subject=$(echo "$alignment" | awk -F "\t" '{print $1}')
	sstart=$(echo "$alignment" | awk -F "\t" '{print $4}')
	send=$(echo "$alignment" | awk -F "\t" '{print $4}')
	info=$(echo "$alignment" | awk -F "\t" '{print $9}')

	query=$(echo "$info" | awk -F ";" '{print $2}' | sed 's/Name=//')
	cov=$(echo "$info" | awk -F ";" '{print $4}' | sed 's/coverage=//')
	pid=$(echo "$info" | awk -F ";" '{print $5}' | sed 's/identity=//')
	matches=$(echo "$info" | awk -F ";" '{print $6}' | sed 's/matches=//')
	mismatches=$(echo "$info" | awk -F ";" '{print $7}' | sed 's/mismatches=//')
	indels=$(echo "$info" | awk -F ";" '{print $8}' | sed 's/indels=//')
	unknowns=$(echo "$info" | awk -F ";" '{print $9}' | sed 's/unknowns=//')

	echo -e "$query\t$subject\t$pid\t$matches\t$mismatches\t$indels\t$unknowns\t$cov\t$sstart\t$send" >> $outfile
	
done

