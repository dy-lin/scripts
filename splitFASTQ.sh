#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $PROGRAM <GZIPPED FASTQ file>" 1>&2
	echo "DESCRIPTION: Takes a gzipped FASTQ file as input, and separates the interleaved reads into two separate files." 1>&2
	exit 1
fi

fastq=$1
filename=$(basename $fastq ".fq.gz")
set -eu -o pipefail

seqtk seq -1 $fastq | sed '/^@/ s/$/\/1/' |pigz -p 64 > ${filename}.R1.fq.gz
seqtk seq -2 $fastq | sed '/^@/ s/$/\/2/' |pigz -p 64 > ${filename}.R2.fq.gz

