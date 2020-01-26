#!/bin/bash

if [[ "$#" -ne 1 ]]
then
	echo "$(basename $0) <CSV or TSV file>" 1>&2
	exit 1
fi

file=$1

ext=${file##*.}

if [[ "$ext" == "tsv" ]]
then
	cat $file | column -t -s "\t" | less -S
fi

if [[ "$ext" == "csv" ]]
then
	cat $file | awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", " ", $i) } 1' | column -t -s "," | less -S
fi
