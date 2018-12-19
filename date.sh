#!/bin/bash

if [[ "$#" -lt 1  || "$#" -gt 2 ]]
then
	echo "USAGE: $0 [DATE] <FILES or PATTERN>"
	echo "Example: $0 DDMMMYYYY *"
	echo "[DATE] is optional-- if one is not provided, today's date is used by default."
	exit 1
fi

# If only one argument is given, use today's date by default
if [ "$#" -eq 1 ]
then
	echo "One argument detected. Using today's date by default."
	date=$(date | awk '{print $3 $2 $6}')
# If two arguments are detected, using the provided date instead
else
	date=$1
	shift
	echo "Two arguments detected. Using $date as date."
fi

# Adding the date to the filename
for file in "$@"
do	
	ext=${file##*.}
	mv "$file" "${file/$ext/$date.$ext}"
done

