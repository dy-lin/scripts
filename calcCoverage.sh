#!/bin/bash
set -eu -o pipefail
if [ "$#" -lt 1 ]
then
	echo "USAGE: $(basename $0) <FASTA file> <reads> <reads>" 1>&2
	exit 1
fi
alignment=$1

# Remove the file extension
temp=$(basename $alignment)
filename=${temp%%.*}

# If the BAM file is already sorted, go ahead
if [[ "$#" -eq 1 && "$alignment" == *sorted*.bam ]]
then
	result=$(median.R <(awk '{print $3}' <(samtools depth -a $alignment)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo -e "Assembly: ${filename}\tCoverage: ${median}x"
# If the BAM file is unsorted, sort it first
elif [[ "$#" -eq 1 && "$alignment" == *.bam ]]
then
	samtools sort $alignment > ${filename}.reads.sorted.bam
	result=$(median.R <(awk '{print $3}' <(samtools depth -a ${filename}.reads.sorted.bam)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo -e "Assembly: ${filename}\tCoverage: ${median}x"

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

	if [[ "$readsf" != *.f?.gz || "$readsb" != *.f?.gz ]]
	then
		echo "Three arguments detected. Invalid reads files." 1>&2
		echo "USAGE: $(basename $0) <FASTA file> <reads> <reads>" 1>&2
		exit 1
	fi
	if [[ ! -e "${assembly}.sa" ]]
	then
		bwa index $assembly 2> /dev/null
	fi
	bwa mem -t 48 $assembly $readsf $readsb 2> /dev/null | samtools view -F 4 -b | samtools sort > ${filename}.reads.sorted.bam
	result=$(median.R <(awk '{print $3}' <(samtools depth -a ${filename}.reads.sorted.bam)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo -e "Assembly: ${filename}\tCoverage: ${median}x"

# If two, assume first is assembly second is interleaved reads.
elif [ "$#" -eq 2 ]
then
	assembly=$1
	reads=$2
	temp=$(basename $assembly)
	filename=${temp%.*}
	if [[ "$assembly" != *.fa && "$assembly" != *.fasta && "$assembly" != *.f?.gz ]]
	then
		echo "Two arguments detected. Invalid assembly file." 1>&2
		echo "USAGE: $(basename $0) <FASTA file> <inter-leaved reads>" 1>&2
		exit 1
	fi
	if [[ "$reads" != *.f?.gz ]]
	then
		echo "Two arguments detected. Invalid reads files." 1>&2
		echo "USAGE: $(basename $0) <FASTA file> <inter-leaved reads>" 1>&2
		exit 1
	fi
	if [[ ! -e "${assembly}.sa" ]]
	then
		bwa index $assembly 2> /dev/null
	fi
	bwa mem -t 48 $assembly $reads 2> /dev/null | samtools view -F 4 -b | samtools sort > ${filename}.reads.sorted.bam
	result=$(median.R <(awk '{print $3}' <(samtools depth -a ${filename}.reads.sorted.bam)) 2> /dev/null)
	median=$(echo $result | awk '{print $2}')
	echo -e "Assembly: ${filename}\tCoverage: ${median}x"

# Or the arguments are invalid
else
	echo "USAGE: $(basename $0) <BAM file (preferably sorted)> OR" 1>&2
	echo "USAGE: $(basename $0) <FASTA file> <reads> <reads>" 1>&2
	echo "DESCRIPTION: Takes a sorted BAM file or a FASTA file with assembly reads and outputs the coverage." 1>&2
	exit 1
fi
