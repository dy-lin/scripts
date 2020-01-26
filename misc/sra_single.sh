#!/bin/bash

set -euo pipefail 
if [[ "$#" -ne 1 ]]
then
	echo "$(basename $0) <SRA accession>" 1>&2
	exit 1
fi

sra=$1
threads=8

fasterq-dump --split-files --threads $threads --progress --force $sra &> ${sra}.log
rm /home/dlin/ncbi/public/sra/${sra}.sra.cache &> /dev/null
pigz -p $threads ${sra}_1.fastq
pigz -p $threads ${sra}_2.fastq
