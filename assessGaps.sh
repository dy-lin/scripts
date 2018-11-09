#!/bin/bash
filename=$(basename $1 ".fa")

seqtk seq $1 > ${filename}.seqtk.fa
if [ "$#" -eq 2 ]
then
	findGaps.py ${filename}.seqtk.fa $2
else
	findGaps.py ${filename}.seqtk.fa
fi

column -t -s$'\t' ${filename}.seqtk.gaps.tsv
