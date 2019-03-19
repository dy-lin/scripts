#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM [protein family] <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) and runs sequences against the Pfam database." 1>&2
	exit 1
fi

if [[ ! -f "$1" ]]
then
	fam=true
	prot=$(echo $1 | sed 's/s$//')
	shift
else
	fam=false
fi

for fasta in "$@"
do	
	if [[ ! -e "$fasta" ]]
	then
		echo "$(basename $fasta) does not exist." 1>&2
		continue
	fi
	echo "Running Pfam on $(basename $fasta)..." 1>&2
	filename=${fasta%.*}
	hmmscan --notextw --noali --tblout ${filename}.tbl -o /dev/null /projects/btl/dlin/datasets/Pfam $fasta
	echo "Pfam domains/families written to $(basename "${filename}.tbl")" 1>&2
	if [[ "$fam" = true ]]
	then
		seqtk subseq $fasta <(awk -v var="$prot" 'IGNORECASE = 1 {if($1==var) print $3}' ${filename}.tbl | sort -u ) > ${filename}.${prot}s.faa
		if [[ ! -s "${filename}.${prot}s.faa" ]]
		then
			echo "No sequences with a Pfam $prot domain/family found." 1>&2
			rm ${filename}.${prot}s.faa
		else
			echo "Sequences with a Pfam $prot domain/family written to $(basename "${filename}.${prot}s.faa")" 1>&2
		fi
	fi
done
