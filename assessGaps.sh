#!/bin/bash
if [ "$#" -eq 0 ]
then
	echo "USAGE: $(basename $0) <FASTA file> [GFF file]"
	exit 1
fi

#Get filename without extension
filename=$(basename $1)
extension="${filename##*.}"
filename=$(basename $1 $extension)


#Run seqtk seq in case FASTA sequence is not all on one line (e.g. Downloaded from NCBI)
seqtk seq $1 > ${filename}seqtk.fa

#If a GFF file is given, findGaps.py will show which gene features contain gaps
if [ "$#" -eq 2 ]
then
	findGaps.py ${filename}seqtk.fa $2
else
	findGaps.py ${filename}seqtk.fa
fi

#Print gap assessment to screen neatly
if [ -e "${filename}.seqtk.gaps.tsv" ]
then
	column -t -s$'\t' ${filename}seqtk.gaps.tsv
fi

#Remove seqtk file as it is no loner needed
rm ${filename}seqtk.fa
