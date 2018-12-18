#!/bin/bash

if [ "$#" -lt 2 ]
then
	echo "USAGE: $0 <DATE> <FILES or PATTERN>"
	echo "Example: $0 DDMMMYYYY *"
	exit 1
fi

date=$1
shift

for file in "$@"
do	
	ext=${file##*.}
	mv "$file" "${file/$ext/$date.$ext}"
done
