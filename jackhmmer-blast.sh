#!/bin/bash
PROGRAM=$(basename $0)
set -eu -o pipefail
if [[ "$#" -ne 3 ]]
then
	echo "USAGE: $PROGRAM <class of AMPs> <database> <output file>" 1>&2
	exit 1
fi

lit=$1
database=$2
outfile=$3

echo "Making BLAST database..."
makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer

echo -e "\nBLASTing..."
blastp -db jackhmmer -query $lit -out jackhmmer.blastp -outfmt '6 std qcovs' -num_threads 48
echo "Running seqtk..."
threshold=90
aligned=$threshold
while [[ ! -s jackhmmer-blast-hits.faa && "$threshold" -ge 50 ]]
do
	seqtk subseq $database <(awk -v var=$threshold '{if ($3>var) print $2}' jackhmmer.blastp | sort -u) > jackhmmer-blast-hits.faa
	aligned=$threshold
	threshold=$((threshold-10))
done

if [[ -s jackhmmer-blast-hits.faa ]]
then
	echo "Proteins aligned with $aligned% identity or higher."
else
	echo "There are no proteins that align with 50% identity or higher."
	echo "Loading jackhmmer-blast-hits.faa with ALL blast hits."
	seqtk subseq $database <(awk '{print $2}' jackhmmer-blastp | sort -u) > jackhmmer-blast-hits.faa
fi
# Clean directories
if [[ "$(basename $(pwd))" != bs* ]]
then
	dir=$(echo $outfile | awk -F "_" '{print $2}')
	mkdir -p $dir

	for i in jackhmmer*
	do
		if [[ -L "$i" ]]
		then
			cp $(readlink -f $i) $dir
		else
			mv $i $dir
		fi
	done
fi
echo "Status: Success."
