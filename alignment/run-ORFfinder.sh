#!/bin/bash
PROGRAM=$(basename $0)
set -euo pipefail
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) and outputs both ORFs and CDS." 1>&2
	exit 1 
fi

mkdir -p log
for fasta in $@
do
	if [[ -e "$fasta" ]]
	then
		if [[ "$(readlink -f $fasta)" == *.gz ]]
		then
			name=$(basename $fasta ".fa.gz")

			if [[ ! -e "${name}.ORFs.faa" ]]
			then
				echo "Running ORFfinder for gzipped $(basename $fasta), outputting an ORF FASTA file..." 1>&2
				ORFfinder -in <(zcat $fasta) -outfmt 0 2>> log/${name}.log | seqtk seq - > ${name}.ORFs.faa
			fi
			if [[ ! -e "${name}.cds.fa" ]]
			then
				echo "Running ORFfinder for gzipped $(basename $fasta), outputting a CDS FASTA file..." 1>&2
				ORFfinder -in <(zcat $fasta) -outfmt 1 2>> log/${name}.log | seqtk seq - > ${name}.cds.fa
			fi
		else
			name=$(basename $fasta ".fa")
			if [[ ! -e "${name}.ORFs.faa" ]]
			then
				echo "Running ORFfinder for $(basename $fasta), outputting an ORF FASTA file..." 1>&2
				ORFfinder -in $fasta -outfmt 0 2>> log/${name}.log | seqtk seq - > ${name}.ORFs.faa
			fi
			if [[ ! -e "${name}.cds.fa" ]]
			then
				echo "Running ORFfinder for $(basename $fasta), outputting a CDS FASTA file..." 1>&2
				ORFfinder -in $fasta -outfmt 1 2>> log/${name}.log | seqtk seq - > ${name}.cds.fa
			fi
		fi
	else
		echo "$(basename $fasta) does not exist." 1>&2
	fi
done



echo "...Done." 1>&2
