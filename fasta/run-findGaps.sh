#!/bin/bash
set -eu -o pipefail
if [ "$#" -eq 0 ]
then
	echo "USAGE: $(basename $0) <FASTA file> [GFF file]" 1>&2
	echo "DESCRIPTION: Takes a FASTA file (and optionally a GFF file) and outputs a gap summary table and a scaftig FASTA file. If a GFF file is also given, features containing gaps are also listed." 1>&2
	exit 1
fi

# Get filename without extension
filename=$(basename $1)
extension="${filename##*.}"
filename=$(basename $1 $extension)

# Run seqtk seq in case FASTA sequence is not all on one line (e.g. Downloaded from NCBI)

if [[ "$(wc -l $1 | awk '{print $1}')" -ne 2 ]]
then
	fileout=${filename}seqtk.fa
	seqtk seq $1 > $fileout
else
	fileout=$1
fi
# If a GFF file is given, findGaps.py will show which gene features contain gaps
if [ "$#" -eq 2 ]
then
	findGaps.py $fileout $2
else
	findGaps.py $fileout
fi

# Print gap assessment to screen neatly
if [ -e "${filename}.seqtk.gaps.tsv" ]
then
	column -t -s$'\t' ${filename}.gaps.tsv
fi

if [[ -e ${filename}seqtk.fa ]]
then
	rm ${filename}seqtk.fa
fi
