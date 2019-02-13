#!/bin/bash
PROGRAM=$(basename $0)
gethelp=false
suffix=false
while getopts :hs opt
do
	case $opt in
		h) gethelp=true;;
		s) suffix=true;;
		\?) echo "$PROGRAM: Invalid option $opt" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))
if [[ "$#" -ne 1 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <GZIPPED FASTQ file>" 1>&2
	echo "DESCRIPTION: Takes a gzipped FASTQ file as input, and separates the interleaved reads into two separate files." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-s\tAdd /1 and /2 suffixes" 1>&2
	exit 1
fi

fastq=$1
filename=$(basename $fastq ".fq.gz")
set -eu -o pipefail

if [[ "$suffix" = true ]]
then 
	seqtk seq -1 $fastq | sed '/^@/ s/$/\/1/' |pigz -p 64 > ${filename}.R1.fq.gz
	seqtk seq -2 $fastq | sed '/^@/ s/$/\/2/' |pigz -p 64 > ${filename}.R2.fq.gz
else
	seqtk seq -1 $fastq | pigz -p 64 > ${filename}.R1.fq.gz
	seqtk seq -2 $fastq | pigz -p 64 > ${filename}.R2.fq.gz
fi


