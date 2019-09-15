#!/bin/bash
set -euo pipefail
# Parses GFF to see old name and new name

gffs="$@"

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <GFF file(s)>" 1>&2
	echo "DESCRIPTION: Takes GFF file(s) as input and outputs the current ID and old alias of each mRNA." 1>&2
	exit 1
fi

for gff in $gffs
do
	echo "$(basename $gff):"
	mrna=$(awk -F "\t" '/\tmRNA\t/ {print $9}' $gff | tr -d ' ')
	for line in $mrna
	do
		id=$(echo "$line" | awk -F "ID=" '{print $2}' | awk -F ";" '{print $1}')
		alias=$(echo "$line" | awk -F "Alias=" '{print $2}' | awk -F ";" '{print $1}')
		echo -e "$id\t$alias"
	done
	echo
done
