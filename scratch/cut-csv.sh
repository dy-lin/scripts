#!/bin/bash

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <CSV or TSV FILE>" 1>&2
	exit 1
fi

CUT=/projects/btl/dlin/bin/backup/miniconda3/bin/csvcut

file=$1
dir=$(dirname $file)

response=""

if [[ "$file" == *.csv ]]
then
	while [[ "$response" == "" ]]
	do
		cat $file | awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", " ", $i) } 1' | column -t -s "," | less -S
		read -p "Please enter the column names you wish to extract (separated by commas), or X to exit: " response
	done
	
	if [[ "$response" == "X" || "$response" == "x" || "$response" == "EXIT" || "$response" == "exit" || "$response" == "Exit" ]]
	then
		exit 1
	else
		${CUT} -c "$response" $file | tr "," "\t" > $dir/treatments.tsv
	fi
elif [[ "$file" == *.tsv ]]
then
	while [[ "$response" == "" ]]
	do
		cat $file | awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", " ", $i) } 1' | column -t -s "\t" | less -S
		read -p "Please enter the column names you wish to extract (separated by commas), or X to exit: " response
	done

	if [[ "$response" == "X" || "$response" == "x" || "$response" == "EXIT" || "$response" == "exit" || "$response" == "Exit"  ]]
	then
		exit 1
	else
		${CUT} -t -c "$response" $file | tr "," "\t" > $dir/treatments.tsv
	fi
else
	echo "Invalid file extension." 1>&2
	exit 1
fi
