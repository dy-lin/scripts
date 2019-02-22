#!/bin/bash

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
