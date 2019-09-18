#!/bin/bash
PROGRAM=$(basename $0)
set -euo pipefail
gethelp=false
putN=0
reverse=false
remove=false
circle=false
while getopts :cdhN:r opt
do
	case $opt in
		c) circle=true; remove=false;;
		d) remove=true; circle=false;;
		h) gethelp=true;;
		N) putN=$OPTARG;;
		r) reverse=true;;
		\?) echo "$PROGRAM: invalid option $OPTARG" 1>&2; exit 1;;
	esac
done

shift $((OPTIND-1))


if [[ "$#" -ne 2 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <segment length> <FASTA file>" 1>&2
	echo "DESCRIPTION: Takes sequence of specified length from the beginning of the sequence and stitches it onto the end of the sequence." 1>&2
	echo -e "OPTIONS:\n\t-c\t\tcircularize by taking <INT> bp from the end and stitching it to <INT> bp from the beginning\n\t-d\t\tdelete the segment\n\t-h\t\tShow help menu\n\t-N <INT>\tadd <INT> number of N's in the stitch\n\t-r\t\treverse (end sequence instead)" 1>&2
	exit 1
fi
begin=1
end=$1
fasta=$2

# Make sure FASTA is a one liner
processFASTA.sh $fasta

length=$(seqtk comp $fasta | awk '{print $2}')
name=$(basename $fasta ".fasta")
name=$(basename $fasta ".fa")
if [[ "$circle" = true ]]
then
	if [[ "$putN" -eq 0 ]]
	then
		echo "ERROR: In order to circularize a sequence, -N <INT> must be used." 1>&2
		echo "USAGE: $PROGRAM <segment length> <FASTA file>" 1>&2
		echo "DESCRIPTION: Takes sequence of specified length from the beginning of the sequence and stitches it onto the end of the sequence." 1>&2
		echo -e "OPTIONS:\n\t-c\t\tcircularize by taking <segment length> bp from the end and stitching it to <segment length> bp from the beginning\n\t-d\t\tdelete the segment\n\t-h\t\tShow help menu\n\t-N <INT>\tadd <INT> number of N's in the stitch\n\t-r\t\treverse (end sequence instead)" 1>&2
		exit 1
	fi
	if [[ "$reverse" = true ]]
	then
		echo "ERROR: -r does not apply when using -c." 1>&2
		echo "UPDATE: Applying -c only." 1>&2
		reverse=false
	fi
	echo "$(head -n1 $fasta) Stitch: $((length-end+1)):$length-----${putN}N-----$begin:$end" > ${name}.circular.fa
	tail -n1 <(run-slice.sh $((length-end+1)) $fasta) | tr -d '\n' >> ${name}.circular.fa
	for i in $(seq 1 $putN)
	do
		echo -n "N" >> ${name}.circular.fa
	done
	tail -n1 <(run-slice.sh $begin $end $fasta) >> ${name}.circular.fa
	exit 0
fi


if [[ "$remove" = true ]]
then
	if [[ "$reverse" = true ]]
	then
		run-slice.sh $begin $((length-end)) $fasta > ${name}.trimmed.fa
	else
		run-slice.sh $((end+1)) $fasta > ${name}.trimmed.fa
	fi
		exit 0
fi

if [[ "$reverse" = true ]]
then
	end=$((length-end))
fi
if [[ "$putN" -gt 0 ]]
then
	echo "Stitching segments together with ${putN}N's..." 1>&2
	# Write a new header
	# Old header with Stitch Coords and Number of N's
	echo "$(head -n1 $fasta) Stitch: $((end+1)):${length}-----${putN}N-----${begin}:${end}" > ${name}.stitched.fa
	tail -n1 <(run-slice.sh $((end+1)) $fasta) | tr -d '\n' >> ${name}.stitched.fa
	for i in $(seq 1 $putN)
	do
		echo -n "N" >> ${name}.stitched.fa
	done
	tail -n1 <(run-slice.sh $begin $end $fasta) >> ${name}.stitched.fa
	echo "End Product: $((end+1)):${length}-----${putN}N-----${begin}:${end}" 1>&2
else
	echo "Stitching segments together..." 1>&2
	echo "$(head -n1 $fasta) Stitch: $((end+1)):${length}----------${begin}:${end}" > ${name}.stitched.fa
	tail -n1 <(run-slice.sh $((end+1)) $fasta) | tr -d '\n' >> ${name}.stitched.fa
	tail -n1 <(run-slice.sh $begin $end $fasta) >> ${name}.stitched.fa
	echo "End Product: $((end+1)):${length}----------${begin}:${end}" 1>&2
fi

