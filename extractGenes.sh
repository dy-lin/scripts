#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
gethelp=false
fileout=false
newname=false
output=genes.fa
while getopts :hfo:n opts
do
	case $opts in
		h) gethelp=true;;
		f) fileout=true;;
		o) output=$OPTARG;fileout=true;;
		n) newname=true;;
		\?) echo "$PROGRAM: invalid option: $OPTARG" >&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 2 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <GFF file> <FASTA file>" >&2
	echo "DESCRIPTION: Takes a GFF file (containing genes from one or more scaffolds), extracts gene coordinates, and then aligns to corresponding scaffold." >&2
	echo -e "OPTIONS:\n\t-h\t\tShow help menu\n\t-f\t\tWrite to file 'genes.fa'\n\t-n\t\tUse the gene name as sequence ID instead of using the scaffold and slice coordinates\n\t-o <FILENAME>\tWrite to file <FILENAME> specified" >&2
	exit 1
fi

gff=$1
scaffolds=$2
IFS=$'\n'
for line in $(awk '/\tgene\t/ {print $1, $4, $5, $7, $9}' $gff)
do
	scaffold=$(echo $line | awk '{print $1}')
	begin=$(echo $line | awk '{print $2}')
	end=$(echo $line | awk '{print $3}')
	strand=$(echo $line | awk '{print $4}')
	gene=$(echo $line | awk '{print $5}' | awk -F ";" '{print $1}' | awk -F "=" '{print $2}')

	if [[ "$strand" == "+" ]]
	then
		# run-slice.sh prints the contents of a FASTA file
		# take tail to only get the sequence
		# write gene name manually
		if [[ "$fileout" = true ]]
		then
			if [[ "$newname" = true ]]
			then
				echo ">${gene}" >> $output
				run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) | tail -n 1 >> $output
			else
				run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) >> $output
			fi
		else
			if [[ "$newname" = true ]]
			then
				echo ">${gene}"
				run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) )| tail -n 1
			else
				run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) 
			fi
		fi
	else
		if [[ "$fileout" = true ]]
		then
			if [[ "$newname" = true ]]
			then
				echo ">${gene}" >> $output
				seqtk seq -r <(run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) ) | tail -n 1 >> $output
			else
				seqtk seq -r <(run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) ) >> $output
			fi
		else
			if [[ "$newname" = true ]]
			then
				echo ">${gene}"
				seqtk seq -r <(run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) ) | tail -n 1
			else
				seqtk seq -r <(run-slice.sh $begin $end <(seqtk subseq $scaffolds <(echo $scaffold) ) )
			fi
		fi
	fi
done
if [[ "$fileout" = true ]]
then
	echo "FASTA sequence written to $output." >&2
fi
