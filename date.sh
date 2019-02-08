#!/bin/bash
set -eu -o pipefail
if [[ "$#" -lt 1  || "$#" -gt 2 ]]
then
	echo "USAGE: $(basename $0) [DATE] <FILES or PATTERN>" 1>&2
	echo "EXAMPLE: $(basename $0) DDMMMYYYY *" 1>&2
	echo "DESCRIPTION: Takes a date and files/pattern, appends that date to those filenames (or the filenames that match the patern). [DATE] is optional-- if one is not provided, today's date is used by default." 1>&2
	exit 1
fi

# If only one argument is given, use today's date by default
if [ "$#" -eq 1 ]
then
	echo "One argument detected. Using today's date by default." 1>&2
	date=$(date | awk '{if($3<10) {print "0" $3 $2 $6} else { print $3 $2 $6}}')
	# print 0 for double digit numeric date
# If two arguments are detected, using the provided date instead
else
	date=$1
	shift
	echo "Two arguments detected. Using $date as date." 1>&2
fi

# Adding the date to the filename
for file in "$@"
do	
	ext=${file##*.}
	mv "$file" "${file/$ext/$date.$ext}"
done

