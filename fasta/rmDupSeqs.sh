#!/bin/bash
set -euo pipefail
PROGRAM=$(basename $0)
if [[ "$#" -lt 2 ]]
then
	echo "USAGE: $PROGRAM <output filename> <FASTA file(s)>" 1>&2
	echo "DESCRIPTION: Takes a FASTA file and removes duplicate sequences, based on sequence ID and the sequence itself. If given more than one FASTA file, concatenates the files and THEN removes duplicate sequences." 1>&2
	exit 1
fi

# First argument must be the output filename
# Make sure that it does not already exist!
outfile=$1
shift
while [[ -e "$outfile" ]]
do
	# Check if outfile has the same name as input arguments
	used=false
	for file in $@
	do
		if [[ "$outfile" == "$file" ]]
		then
			used=true
			break
		fi
	done
	if [[ "$used" = true ]]
	then
		echo "Output file is one of the input files." 1>&2
		read -p "Enter new output filename: " outfile
	else
		echo "Output file already exists." 1>&2
		read -p "Overwrite it? " response
		case $response in
			[Yy]*) break
				;;
			[Nn]*) read -p "Enter new output filename: " outfile
				;;
		esac
	fi
done

# Feed Through Seqtk seq to ensure one line sequence
# Translate new lines to tabs
# Replace Tabs and > with Newline and >, that way each header and sequence is on a separate line
# sort and unique (no fields specified, use all)
# translate the tabs back to new lines
# remove any empty lines

cat $@ | seqtk seq | tr '\n' '\t' | sed 's/\t>/\n>/g' | sort -u | tr '\t' '\n' | sed '/^$/d' > $outfile

if [[ "$#" -eq 1 ]]
then
	read -p "Overwrite $1 with $outfile? " response
	case $response in
		[Yy]*) mv $outfile $1
			;;
	esac
fi
