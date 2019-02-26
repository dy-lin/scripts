#!/bin/bash
set -eu -o pipefail

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <bit score threshold>" 1>&2
	exit 1
fi

threshold=$1
# Get all AMPs of a class
amp=$(basename $(pwd))
outfile=${amp}_summary.jira
tissues=$(ls -l | awk '/^d/ {print $9}' | sort -u)

max=$(($(echo $tissues | wc -w) +1 ))
echo -n "||$amp||" > $outfile
count=1

# On the last tissue, write with newline.
for tissue in $tissues
do
	count=$((count+1))
	if [[ "$count" -lt "$max" ]]
	then
		echo -n "$tissue||" >> $outfile
	else
		echo "$tissue||" >> $outfile
	fi
done
for prot in $(for file in $(ls */bs${threshold}/transcripts/*); do basename $file; done | sort -u)
do
	echo -n "|$(basename $prot ".transcript.fa")|" >> $outfile
	count=1
	for tissue in $tissues
	do
		count=$((count+1))
		if [[ -s "${tissue}/bs${threshold}/transcripts/${prot}" ]]
		then
			if [[ "$count" -lt "$max" ]]
			then
				echo -n "(/)|" >> $outfile
			else
				echo "(/)|" >> $outfile
			fi
		else
			if [[ "$count" -lt "$max" ]]
			then
				echo -n "(x)|" >> $outfile
			else
				echo "(x)|" >> $outfile
			fi
		fi
	done
done
echo -n "|Total|" >> $outfile
count=1
for tissue in $tissues
do
	num_transcripts=$(ls ${tissue}/bs${threshold}/transcripts/* | wc -l)
	if [[ "$count" -lt "$max" ]]
	then
		echo -n "$num_transcripts|" >> $outfile
	else
		echo "$num_transcripts|" >> $outfile
	fi
done

echo >> $outfile

echo "||tissue||\# transcripts||\# ORFs||jackhmmer hits||blast hits|| corresponding \# transcripts || defensins (including partial) || defensin-like || threshold ||" >> $outfile
for tissue in $tissues
do
	num_transcriptome=$(zgrep -c '^>' ../assemblies/${tissue})
	num_ORFs=$(grep -c '^>' ${tissue}/${tissue}.ORFs.fa)
	num_hits=$(grep -c '^>' ${tissue}/bs${threshold}/jackhmmer-hits.faa)
	num_blast=$(grep -c '^>' ${tissue}/bs${threshold}/jackhmmer-blast-hits.faa)
	num_transcripts=$(ls ${tissue}/bs${threshold}/transcripts/* | wc -l)
	if [[ -e jackhmmer-blast-hits.trimmed.faa ]]
	then
		num_amplike=$(grep -c 'defensin-like' ${tissue}/bs${threshold}/jackhmmer-blast-hits.trimmed.fa)
	else
		num_amplike=0
	fi
	num_amp=$((num_blast-num_amplike))

	echo "|$tissue|$num_transcriptome|$num_ORFs|$num_hits|$num_blast|$num_transcripts|$num_amp|$num_amplike|" >> $outfile
done
