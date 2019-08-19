#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $PROGRAM [prefix] <FASTA files>" 1>&2
	echo "DESCRIPTION: Rename FASTA file(s) with prefix." 1>&2
	echo "EXAMPLE: $PROGRAM assembly.kX-kcX-XXM.fa" 1>&2
	echo -e "\tDefault prefix is XXM-kcX-kX, where header becomes >XXM-kcX-kX-header" 1>&2
	exit 1
fi

if [[ "$#" -eq 1 && "$1" != *.fa && "$!" != *.fasta ]]
then
	echo "No FASTA file provided." 1>&2
	echo "USAGE: $PROGRAM [prefix] <FASTA files>" 1>&2
	echo "DESCRIPTION: Rename FASTA file(s) with prefix. If not prefix is given, the last 3 words before the file extension is used." 1>&2
	echo "EXAMPLE: $PROGRAM assembly.kX-kcX-XXM.fa" 1>&2
	echo -e "\tDefault prefix is XXM-kcX-kX, where header becomes >XXM-kcX-kX-header" 1>&2
	exit 1
fi
if [[ "$1" != *.fa && "$1" != *.fasta ]]
then
	prefix=$1
	shift
	for fasta in $@
	do
		filename=$(echo $fasta | sed 's/.fasta$//' | sed 's/.fa$//')
		sed "s/^>/>$prefix-/" $fasta > ${filename}.prefix.fa
	done
else
	for fasta in $@
	do
		filename=$(echo $fasta | sed 's/.fasta$//' | sed 's/.fa$//')
		prefix_old=$(echo $fasta | awk -F "." '{print $(NF-1)}')
		k=$(echo $prefix_old | awk -F "-" '{print $1}')
		kc=$(echo $prefix_old | awk -F "-" '{print $2}')
		read_pairs=$(echo $prefix_old |awk -F "-" '{print $3}')

		if [[ -z "$k" || -z "$kc" || -z "$read_pairs" ]]
		then
			echo "File is not in XXM-kcX-kX.fa format." 1>&2
			exit 1
		fi
		prefix="$read_pairs-$kc-$k"
		sed "s/^>/>$prefix-/" $fasta > ${filename}.prefix.fa
	done
fi


		

