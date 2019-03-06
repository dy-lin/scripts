#!/bin/bash
PROGRAM=$(basename $0)
lit=$1
database=$2
outfile=$3

if [[ "$#" -ne 3 ]]
then
	echo "USAGE: $PROGRAM <class of AMPs> <database> <output file>" 1>&2
	exit 1
fi

echo "Making BLAST database..."
makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer

echo -e "\nBLASTing..."
blastp -db jackhmmer -query $lit -out jackhmmer.blastp -outfmt '6 std qcovs' -num_threads 48
echo "Running seqtk..."
seqtk subseq $database <(awk '{if ($3>90) print $2}' jackhmmer.blastp | sort -u) > jackhmmer-blast-hits.faa

# Clean directories
dir=$(echo $outfile | awk -F "_" '{print $2}')
mkdir -p $dir
mv jackhmmer* $dir
echo "Status: Success."
