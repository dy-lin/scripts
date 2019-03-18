#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <protein family> <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) and runs sequences against the Pfam database." 1>&2
	exit 1
fi
prot=$(echo $1 | sed 's/s$//')
shift
for fasta in "$@"
do
	filename=${fasta%.*}
	hmmscan --notextw --noali --tblout ${filename}.tbl -o /dev/null /projects/btl/dlin/datasets/Pfam $fasta
	seqtk subseq $fasta <(awk -v var="$prot" 'IGNORECASE = 1 {if($1==var) print $3}' ${filename}.tbl | sort -u ) > ${filename}.${prot}s.faa
done
