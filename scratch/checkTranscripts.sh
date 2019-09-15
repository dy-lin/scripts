#!/bin/bash
set -euo pipefail
if [[ ! -s jackhmmer-blast-hits.faa ]]
then
	echo "USAGE: $(basename $0)" 1>&2
	echo "DESCRIPTION: Checks whether or not each ORF hit is unique to one transcript." 1>&2
	exit 1
fi
hits=$(awk '/^>/ {print $1}' jackhmmer-blast-hits.faa | wc -l)
transcripts=$(awk '/^>/ {print $1}' jackhmmer-blast-hits.faa | awk -F ":" '{print $1}' | awk -F "_" '{print $2}' | sort -u | wc -l)

if [[ "$hits" -ne "$transcripts" ]]
then
	duplicates=$(awk '/^>/ {print $1}' jackhmmer-blast-hits.faa | awk -F ":" '{print $1}' | awk -F "_" '{print $2}' | sort | uniq -c | awk '{if($1>1) print $2}')
	echo "$hits hits originated from $transcripts transcripts."
	dupcount=$(awk '/^>/ {print $1}' jackhmmer-blast-hits.faa | awk -F ":" '{print $1}' | awk -F "_" '{print $2}' | sort | uniq -c | awk '{if($1>1) print $2}' | wc -l)
	echo "$dupcount transcripts resulted in more than one hit."
	echo
	for t in $duplicates
	do
		echo "Transcript: $t"
		grep -A1 "$t" jackhmmer-blast-hits.faa
		echo
	done
else
	echo "$hits hits originated from $transcripts transcripts."
fi
