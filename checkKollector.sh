#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)

if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $PROGRAM <FILES>" 1>&2
	echo "DESCRIPTION: Takes files as input and checks if the file exists and has a size greater than zero. Intended for checking if Kollector has successfully assembled its targets." 1>&2
	exit 1
fi

# Sort by Run number
files=$(sort -t. -k3 <(echo "$@" | tr ' ' '\n'))

for i in $files
do
	if [[ -s "$i" ]]
	then
		echo "$i exists and has a file size larger than zero." 1>&2
	else
		if [[ -e "$i" ]]
		then
			echo "$i is empty." 1>&2
		else
			echo "$i does not exist." 1>&2
		fi
	fi
done
