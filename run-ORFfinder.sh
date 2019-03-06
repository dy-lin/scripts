#!/bin/bash
PROGRAM=$(basename $0)

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <FASTA file>" 1>&2
	echo "DESCRIPTION: Takes a FASTA file and outputs both ORFs and CDS." 1>&2
	exit 1 
fi

mkdir -p log
for fasta in $@
do
	if [[ "$fasta" == *.gz ]]
	then
		name=$(basename $fasta ".fa.gz")
		ORFfinder -in <(zcat $fasta) -outfmt 0 | seqtk seq - > ${name}.ORFs.faa 2>> log/${name}.log
		ORFfinder -in <(zcat $fasta) -outfmt 1 | seqtk seq - > ${name}.cds.fa 2>> log/${name}.log
	else
		name=$(basename $fasta ".fa")
		ORFfinder -in $fasta -outfmt 0 | seqtk seq - > ${name}.ORFs.faa 2>> log/${name}.log
		ORFfinder -in $fasta -outfmt 1 | seqtk seq - > ${name}.cds.fa 2>> log/${name}.log
	fi
done




