#!/bin/bash
set -euo pipefail
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <MAP> <GFF>" 1>&2
	exit 1
fi

map=$1
shift

for gff in $@
do
	genes=$(awk -F "\t" '/\tgene\t/ {print $9}' $gff | awk -F ";" '{print $1}' | sed 's/ID=//')
	mrnas=$(awk -F "\t" '/\tmRNA\t/ {print $9}' $gff | awk -F ";" '{print $1}' | sed 's/ID=//')
	for mrna in $mrnas
	do
		sub=$(awk -v var="$mrna" '{if($1==var) print $2}' $map)
		if [[ ! -z "$sub" ]]
		then
			sed -i "s/$mrna/$sub/g" $gff
		fi
	done
	for gene in $genes
	do
		sub=$(awk -v var="$gene" '{if($1==var) print $2}' $map)
		if [[ ! -z "$sub" ]]
		then
			sed -i "s/$gene/$sub/g" $gff
		fi
	done
done
