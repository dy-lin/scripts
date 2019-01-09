#!/bin/bash

gff=$1

# Split a merged GFF of various scaffolds (e.g. of one spruce genotype) into its respective GFF files for individual viewing in IGV

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <GFF file>" 1>&2 
	echo "DESCRIPTION: Takes a GFF file of multiple scaffolds and splits them by scaffold." 1>&2
	exit 1
fi

csplit -s $gff /##gff-version/ {*}

for file in xx*
do
	if [ "$file" != "xx00" ]
	then
		filename=$(awk 'NR==2' $file | awk '{print $2}')
		mv $file ${filename}.gff
	fi
done

rm "xx00"
