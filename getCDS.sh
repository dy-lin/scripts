#!/bin/bash

scaf=$1
gff=$2
name=$(basename "${gff%%.*}")
# Fetches the CDS separately for each GFF entry

if [[ -e "${name}.temp.fa" ]]
then
	rm "${name}.temp.fa"
fi

gt extractfeat -retainids -matchdescstart -type CDS -seqfile $scaf -o ${name}.temp.fa $gff

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
	grep -A1 "$header" ${name}.temp.fa | grep -v '^>' | tr -d '\n' >> ${name}.cds.fa
	echo >> ${name}.cds.fa
done

rm ${name}.temp.fa
