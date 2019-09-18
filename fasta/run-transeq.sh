#!/bin/bash
set -euo pipefail
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <nucleotide sequence(s) or FASTA file(s)>" 1>&2
	exit 1
fi
infile=false
if [[ "$infile" = false ]]
then
	# check that arguments are not files
	allFiles=true
	for file in $@
	do
		if [[ ! -f "$file" ]]
		then
			allFiles=false
		fi
	done

	if [[ "$allFiles" = true ]]
	then
		infile=true
	fi
fi


if [[ "$infile" = false ]]
then
	# if no infile, make a temp file to work with.
	count=1
	filenum=1
	outfile=temp.fa
	while [[ -e "$outfile" ]]
	do
		filenum=$((filenum+1))
		outfile="temp${filenum}.fa"
	done

	for seq in $@
	do
		if [[ ! -f "$seq" ]]
		then
			echo -e ">header$count\n$seq" >> $outfile
			count=$((count+1))
		else
			cat $seq >> $outfile
		fi
	done
else
	if [[ "$#" -gt 1 ]]
	then
		count=1
		filenum=1
		outfile=temp.fa
		while [[ -e "$outfile" ]]
		do
			filenum=$((filenum+1))
			outfile="temp${filenum}.fa"
		done
		
		for file in $@
		do
			if [[ -e "$file" ]]
			then
				cat $file >> $outfile
			else
				echo "$(basename $file) does not exist." 1>&2
				rm $outfile
				exit 1
			fi
		done

	else
		if [[ -e "$1" ]]
		then
			outfile=$1
		else
			echo "$(basename $1) does not exist." 1>&2
			exit 1
		fi
	fi
fi
transeq -sformat pearson -sequence $outfile -outseq /dev/stdout 2> /dev/null | sed 's/_1$//' | seqtk seq -
if [[  "$outfile" == temp* ]]
then
	rm $outfile
fi

