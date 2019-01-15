#!/bin/bash
set -eu -o pipefail
alignment=$1

# Remove the file extension
temp=$(basename $alignment)
filename=${temp%%.*}

# If the BAM file is already sorted, go ahead
if [[ "$#" -eq 1 && "$alignment" == *sorted*.bam ]]
then
	result=$(median.R <(awk '{print $3}' <(samtools depth -a $alignment)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo "Coverage: ${median}x"
# If the BAM file is unsorted, sort it first
elif [[ "$#" -eq 1 && "$alignment" == *.bam ]]
then
	samtools sort $alignment > ${filename}.sorted.bam
	result=$(median.R <(awk '{print $3}' <(samtools depth -a ${filename}.sorted.reads.bam)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo "Coverage: ${median}x"
# If there are three arguments, assume they are assembly and reads
elif [ "$#" -eq 3 ]
then
	assembly=$1
	readsf=$2
	readsb=$3
	temp=$(basename $assembly)
	filename=${temp%.*}
	# Check if arguments have the correct file extensions
	if [[ "$assembly" != *.fa && "$assembly" != *.fasta && "$assembly" != *.f?.gz ]]
	then
		echo "Three arguments detected. Invalid assembly file." 1>&2
		echo "USAGE: $(basename $0) <FASTA file> <reads> <reads>" 1>&2
		exit 1
	fi

	if [[ "$readsf" != *R?.f?.gz || "$readsb" != *R?.f?.gz ]]
	then
		echo "Three arguments detected. Invalid reads files." 1>&2
		echo "USAGE: $(basename $0) <FASTA file> <reads> <reads>" 1>&2
		exit 1
	fi
	bwa index $assembly
	bwa mem -t 16 $assembly $readsf $readsb | samtools view -F 4 -b -o ${filename}.reads.bam
	samtools sort ${filename}.reads.bam > ${filename}.reads.sorted.bam
	result=$(median.R <(awk '{print $3}' <(samtools depth -a ${filename}.reads.sorted.bam)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo "Coverage: ${median}x"
# Or the arguments are invalid
else
	echo "USAGE: $(basename $0) <BAM file (preferably sorted)> OR" 1>&2
	echo "USAGE: $(basename $0) <FASTA file> <reads> <reads>" 1>&2
	echo "DESCRIPTION: Takes a sorted BAM file or a FASTA file with assembly reads and outputs the coverage." 1>&2
	exit 1
fi
