#!/bin/bash

# Run volcanoDE.R for the outputfile.
if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <DE gene file> <InterProScan TSV file>" 1>&2
	echo "DESCRIPTION: Creates a reference GOterm file and extracts the significant DE genes." 1>&2
	exit 1
fi

infile=$1
ips=$2


kallisto_dir=$(echo "$infile" | sed 's/kallisto.\+$/kallisto/')
specific_dir=$(dirname $infile)
filename=$(basename $infile ".txt")

awk -F "\t" 'BEGIN{OFS="\t"}{print $1, $14}' $ips | grep 'GO:' > $kallisto_dir/gene-GOterms.tsv
grep -f $infile $kallisto_dir/gene-GOterms.tsv > $specific_dir/${filename}.GOterms.tsv


