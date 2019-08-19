#!/bin/bash
PROGRAM=$(basename $0)
set -eu -o pipefail
gethelp=false
verbose=false
while getopts :hv opt
do
	case $opt in
		h) gethelp=true;;
		v) verbose=true;;
		\?) echo "ERROR: $PROGRAM: Invalid option $opt" 1>&2 ; exit 1;;
	esac
done
command="COMMAND: $PROGRAM $*"
shift $((OPTIND-1))
if [[ "$#" -ne 3  ]]
then
	echo "USAGE: $PROGRAM <AMP sequences> <JACKHMMER database> <JACKHMMER output file>" 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-v\tverbose logging" 1>&2
	exit 1
fi

lit=$1
database=$2
outfile=$3
if [[ "$verbose" = true ]]
then
	echo -e"\t$command" 1>&2
fi
echo "Making BLAST database..." 1>&2
if [[ "$verbose" = true ]]
then
	echo -e "\tCOMMAND:  makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer" 1>&2
fi
makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer

echo -e "\nBLASTing..."
if [[ "$verbose" = true ]]
then
	echo -e "\tCOMMAND: blastp -db jackhmmer -query $lit -out jackhmmer.blastp -outfmt '6 std qcovs' -num_threads 48" 1>&2
fi

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
