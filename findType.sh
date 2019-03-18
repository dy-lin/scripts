#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) and outputs 'nucl' or 'prot', also renaming the file extension." 1>&2
	exit 1
fi

for fasta in "$@"
do
	seq=$(awk '!/^>/' $fasta | tr -d '\n')
	res=$(residue.py "$seq")
	case $res in
		prot) filename=${fasta%.*};mv $fasta ${filename}.faa;echo "$res";;
		nucl) echo "$res";;
		invalid) echo "$res"; echo "Invalid residue present." 1>&2;;
		ambiguous) echo "$res"; echo "Ambiguous base present." 1>&2 ;;
	esac
done



