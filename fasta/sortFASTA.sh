#!/bin/bash
PROGRAM=$(basename $0)
set -euo pipefail
inplace=false
gethelp=false

while getopts :hi opt
do
	case $opt in
		i) inplace=true;;
		h) gethelp=true;;
		\?) echo "$PROGRAM: invalid option: $opt" 1>&2; exit 1;;
	esac
done
shift "$((OPTIND-1))"

if [[ "$#" -lt 1 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <FASTA file>" 1>&2
	echo "DESCRIPTION: Takes FASTA file(s) as input, and sorts each one alphabetically by ID." 1>&2
	echo "OPTIONS:\n\t-i\tSort the FASTA file(s) in place\n\t-h\tShow help menu" 1>&2
	exit 1
fi


count=0
IFS=$'\n'
for file in "$@"
do
	count=$((count+1))
	if [[ "$file" == *.gz ]]
	then
		echo "FASTA file compressed with gzip detected..." 1>&2
		if [[ "$file" == *.fasta.gz ]]
		then
			filename=$(basename "$file" ".fasta.gz")
			extension=".fasta"
		elif [[ "$file" == *.fa.gz ]]
		then
			filename=$(basename "$file" ".fa.gz")
			extension=".fa"
		else
			echo "ERROR: $(basename $file) is not a FASTA file." 1>&2
			if [[ "$count" -ne "$#" ]]
			then
				echo "Continuing with other files..." 1>&2
				continue
			else
				echo "Try again." 1>&2
				exit 1
			fi
		fi
		if [[ -e ${filename}.sorted${extension} ]]
		then
			rm ${filename}.sorted${extension}
		fi
		for sequence in $(awk '/^>/' <(zcat $file) | sort)
		do
			grep -A1 --no-group-separator "$sequence" $file >> ${filename}.sorted${extension}
		done
		if [[ "$inplace" = true ]]
		then
			mv ${filename}.sorted${extension} $file
			gzip $file
		else
			gzip ${filename}.sorted${extension}
		fi
	else
		if [[ "$file" == *.fasta ]]
		then
			filename=$(basename "$file" ".fasta")
			extension=".fasta"
		elif [[ "$file" == *.fa ]]
		then
			filename=$(basename "$file" ".fa")
			extension=".fa"
		else
			echo "ERROR: $(basename $file) is not a FASTA file." 1>&2
			if [[ "$count" -ne "$#" ]]
			then
				echo "Continuing with other files..." 1>&2
				continue
			else
				echo "Try again." 1>&2
				exit 1
			fi
		fi
		if [[ -e ${filename}.sorted${extension} ]]
		then
			rm ${filename}.sorted${extension}
		fi
		for sequence in $(awk '/^>/' $file | sort)
		do
			grep -A1 --no-group-separator "$sequence" $file >> ${filename}.sorted${extension}
		done
		if [[ "$inplace" = true ]]
		then
			mv ${filename}.sorted${extension} $file
		fi
	fi
done
