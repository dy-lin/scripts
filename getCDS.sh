#!/bin/bash
if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <scaffold> <GFF file>" 1>&2
	echo "DESCRIPTION: Extracts the CDS sequence from the scaffold and GFF annotation file." 1>&2
	exit 1 
fi
scaf=$1
gff=$2
name=$(basename "${gff%%.*}")
# Fetches the CDS separately for each GFF entry

if [[ -e "${name}.temp.fa" ]]
then
	rm "${name}.temp.fa"
fi

gt extractfeat -retainids -matchdescstart -type CDS -seqfile $scaf -o ${name}.temp.fa $gff 2> /dev/null
if [[ "$?" -eq 1 ]]
then
	gt gff3 -sort -retainids $gff > ${name}.sorted.gff
	if [[ "$?" -eq 0 ]]
	then
		gt extractfeat -retainids -matchdescstart -type CDS -seqfile $scaf -force -o ${name}.temp.fa ${name}.sorted.gff
	else
		exit 1
	fi
fi

# Get unique FASTA headers
# For each one, seqtk subseq, then combine together to write a new one

# If ${name}.cds.fa already exists, delete
if [[ -e "${name}.cds.fa" ]]
then
	rm "${name}.cds.fa"
fi

for header in $(grep '^>' ${name}.temp.fa | sort -u )
do
	echo "$header" >> ${name}.cds.fa
#	seqtk subseq ${name}.temp.fa <(echo "$header") | grep -v '^>' | tr -d '\n' >> ${name}.cds.fa
	grep -A1 "$header" ${name}.temp.fa | grep -v '^>' | tr -d '\n' >> ${name}.cds.fa
	echo >> ${name}.cds.fa
done

rm ${name}.temp.fa
