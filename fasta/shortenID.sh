#!/bin/bash
PROGRAM=$(basename $0)
set -euo pipefail
if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) and shortens the IDs to the first word." 1>&2
	exit 1
fi

for fasta in "$@"
do
	filename=$(basename $fasta ".fa")
	filename=$(basename $filename ".fasta")
	while read id
	do
		echo $id | awk '{print $1}' >> ${filename}.shortened.fa
		read seq
		echo $seq >> ${filename}.shortened.fa
	done < $fasta
done
