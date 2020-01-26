#!/bin/bash

if [[ "$#" -eq 0 ]]
then
	echo "$(basename $0) <CSV file(s)>" 1>&2
	exit 1 
fi

file=$1
filename=${file%.*}

cat $file | awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", " ", $i) } 1' > ${filename}.processed.csv
