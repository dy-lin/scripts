#!/bin/bash

# mlr.sh KEYWORD stuff to join

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <keyword> <TSV files>" 1>&2
	exit 1
fi
keyword=$1
shift
files=$(tac <(echo $* | tr ' ' '\n' | sort ))
file2=$(echo "$files" | head -n1)
file1=$(echo "$files" | awk 'NR==2')
base_command="mlr --tsv join -f $file1 -j $keyword $file2"
current_command=$base_command
if [[ "$#" -gt 2 ]]
then
	for file in $(echo "$files" | tail -n +3)
	do
		repeat_command="| mlr --tsv join -f $file -j $keyword"
		current_command="$current_command $repeat_command"
	done
fi

eval $current_command
